/***************************************************************************
For the Third Data Prefetching Championship - DPC3

Paper ID: #4
Instruction Pointer Classifying Prefetcher - IPCP

Authors:
Samuel Pakalapati - samuelpakalapati@gmail.com
Biswabandan Panda - biswap@cse.iitk.ac.in
***************************************************************************/

// Code Refactored by Munawira Kotyad to work with the latest Champsim version

#include "cache.h"

#define NUM_IP_TABLE_L1_ENTRIES 1024 // IP table entries
#define NUM_GHB_ENTRIES 16           // Entries in the GHB
#define NUM_IP_INDEX_BITS 10         // Bits to index into the IP table
#define NUM_IP_TAG_BITS 6            // Tag bits per IP table entry
#define S_TYPE 1                     // stream
#define CS_TYPE 2                    // constant stride
#define CPLX_TYPE 3                  // complex stride
#define NL_TYPE 4                    // next line

// #define SIG_DEBUG_PRINT
#ifdef SIG_DEBUG_PRINT
#define SIG_DP(x) x
#else
#define SIG_DP(x)
#endif

class IP_TABLE_L1
{
public:
  uint64_t ip_tag;
  uint64_t last_page;      // last page seen by IP
  uint64_t last_cl_offset; // last cl offset in the 4KB page
  int64_t last_stride;     // last delta observed
  uint16_t ip_valid;       // Valid IP or not
  int conf;                // CS conf
  uint16_t signature;      // CPLX signature
  uint16_t str_dir;        // stream direction
  uint16_t str_valid;      // stream valid
  uint16_t str_strength;   // stream strength

  IP_TABLE_L1()
  {
    ip_tag = 0;
    last_page = 0;
    last_cl_offset = 0;
    last_stride = 0;
    ip_valid = 0;
    signature = 0;
    conf = 0;
    str_dir = 0;
    str_valid = 0;
    str_strength = 0;
  };
};

class DELTA_PRED_TABLE
{
public:
  int delta;
  int conf;

  DELTA_PRED_TABLE()
  {
    delta = 0;
    conf = 0;
  };
};

IP_TABLE_L1 trackers_l1[NUM_CPUS][NUM_IP_TABLE_L1_ENTRIES];
DELTA_PRED_TABLE DPT_l1[NUM_CPUS][4096];
uint64_t ghb_l1[NUM_CPUS][NUM_GHB_ENTRIES];
uint64_t prev_cpu_cycle[NUM_CPUS];
uint64_t num_misses[NUM_CPUS];
float mpkc[NUM_CPUS] = {0};
int spec_nl[NUM_CPUS] = {0};

/***************Updating the signature*************************************/
uint16_t update_sig_l1(uint16_t old_sig, int delta)
{
  uint16_t new_sig = 0;
  int sig_delta = 0;

  // 7-bit sign magnitude form, since we need to track deltas from +63 to -63
  sig_delta = (delta < 0) ? (((-1) * delta) + (1 << 6)) : delta;
  new_sig = ((old_sig << 1) ^ sig_delta) & 0xFFF; // 12-bit signature

  return new_sig;
}

/****************Encoding the metadata***********************************/
uint32_t encode_metadata(int stride, uint16_t type, int spec_nl)
{

  uint32_t metadata = 0;

  // first encode stride in the last 8 bits of the metadata
  if (stride > 0)
    metadata = stride;
  else
    metadata = ((-1 * stride) | 0b1000000);

  // encode the type of IP in the next 4 bits
  metadata = metadata | (type << 8);

  // encode the speculative NL bit in the next 1 bit
  metadata = metadata | (spec_nl << 12);

  return metadata;
}

/*********************Checking for a global stream (GS class)***************/

void check_for_stream_l1(int index, uint64_t cl_addr, uint8_t cpu)
{
  int pos_count = 0, neg_count = 0, count = 0;
  uint64_t check_addr = cl_addr;

  // check for +ve stream
  for (int i = 0; i < NUM_GHB_ENTRIES; i++) {
    check_addr--;
    for (int j = 0; j < NUM_GHB_ENTRIES; j++)
      if (check_addr == ghb_l1[cpu][j]) {
        pos_count++;
        break;
      }
  }

  check_addr = cl_addr;
  // check for -ve stream
  for (int i = 0; i < NUM_GHB_ENTRIES; i++) {
    check_addr++;
    for (int j = 0; j < NUM_GHB_ENTRIES; j++)
      if (check_addr == ghb_l1[cpu][j]) {
        neg_count++;
        break;
      }
  }

  if (pos_count > neg_count) { // stream direction is +ve
    trackers_l1[cpu][index].str_dir = 1;
    count = pos_count;
  } else { // stream direction is -ve
    trackers_l1[cpu][index].str_dir = 0;
    count = neg_count;
  }

  if (count > NUM_GHB_ENTRIES / 2) { // stream is detected
    trackers_l1[cpu][index].str_valid = 1;
    if (count >= (NUM_GHB_ENTRIES * 3) / 4) // stream is classified as strong if more than 3/4th entries belong to stream
      trackers_l1[cpu][index].str_strength = 1;
  } else {
    if (trackers_l1[cpu][index].str_strength == 0) // if identified as weak stream, we need to reset
      trackers_l1[cpu][index].str_valid = 0;
  }
}

/**************************Updating confidence for the CS class****************/
int update_conf(int stride, int pred_stride, int conf)
{
  if (stride == pred_stride) { // use 2-bit saturating counter for confidence
    conf++;
    if (conf > 3)
      conf = 3;
  } else {
    conf--;
    if (conf < 0)
      conf = 0;
  }

  return conf;
}

void CACHE::prefetcher_initialize() { cout << "CPU " << cpu << " L1D IPCP prefetcher" << endl; }

uint32_t CACHE::prefetcher_cache_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in)
{

  uint64_t curr_page = addr >> LOG2_PAGE_SIZE;
  uint64_t cl_addr = addr >> LOG2_BLOCK_SIZE;
  uint64_t cl_offset = (addr >> LOG2_BLOCK_SIZE) & 0x3F;
  uint16_t signature = 0, last_signature = 0;
  int prefetch_degree = 0;
  int spec_nl_threshold = 0;
  int num_prefs = 0;
  uint32_t metadata = 0;
  uint16_t ip_tag = (ip >> NUM_IP_INDEX_BITS) & ((1 << NUM_IP_TAG_BITS) - 1);

  if (NUM_CPUS == 1) {
    prefetch_degree = 3;
    spec_nl_threshold = 15;
  } else { // tightening the degree and MPKC constraints for multi-core
    prefetch_degree = 2;
    spec_nl_threshold = 5;
  }

  // update miss counter
  if (cache_hit == 0)
    num_misses[cpu] += 1;

  // update spec nl bit when num misses crosses certain threshold
  if (num_misses[cpu] == 256) {
    mpkc[cpu] = ((float)num_misses[cpu] / (current_core_cycle[cpu] - prev_cpu_cycle[cpu])) * 1000;
    prev_cpu_cycle[cpu] = current_core_cycle[cpu];
    if (mpkc[cpu] > spec_nl_threshold)
      spec_nl[cpu] = 0;
    else
      spec_nl[cpu] = 1;
    num_misses[cpu] = 0;
  }

  // calculate the index bit
  int index = ip & ((1 << NUM_IP_INDEX_BITS) - 1);
  if (trackers_l1[cpu][index].ip_tag != ip_tag) { // new/conflict IP
    if (trackers_l1[cpu][index].ip_valid == 0) {  // if valid bit is zero, update with latest IP info
      trackers_l1[cpu][index].ip_tag = ip_tag;
      trackers_l1[cpu][index].last_page = curr_page;
      trackers_l1[cpu][index].last_cl_offset = cl_offset;
      trackers_l1[cpu][index].last_stride = 0;
      trackers_l1[cpu][index].signature = 0;
      trackers_l1[cpu][index].conf = 0;
      trackers_l1[cpu][index].str_valid = 0;
      trackers_l1[cpu][index].str_strength = 0;
      trackers_l1[cpu][index].str_dir = 0;
      trackers_l1[cpu][index].ip_valid = 1;
    } else { // otherwise, reset valid bit and leave the previous IP as it is
      trackers_l1[cpu][index].ip_valid = 0;
    }

    // issue a next line prefetch upon encountering new IP
    uint64_t pf_address = ((addr >> LOG2_BLOCK_SIZE) + 1) << LOG2_BLOCK_SIZE; // BASE NL=1, changing it to 3
    metadata = encode_metadata(1, NL_TYPE, spec_nl[cpu]);
    prefetch_line(ip, addr, pf_address, FILL_L1, metadata);
    return metadata_in;
  } else { // if same IP encountered, set valid bit
    trackers_l1[cpu][index].ip_valid = 1;
  }

  // calculate the stride between the current address and the last address
  int64_t stride = 0;
  if (cl_offset > trackers_l1[cpu][index].last_cl_offset)
    stride = cl_offset - trackers_l1[cpu][index].last_cl_offset;
  else {
    stride = trackers_l1[cpu][index].last_cl_offset - cl_offset;
    stride *= -1;
  }

  // don't do anything if same address is seen twice in a row
  if (stride == 0)
    return metadata_in;

  // page boundary learning
  if (curr_page != trackers_l1[cpu][index].last_page) {
    if (stride < 0)
      stride += 64;
    else
      stride -= 64;
  }

  // update constant stride(CS) confidence
  trackers_l1[cpu][index].conf = update_conf(stride, trackers_l1[cpu][index].last_stride, trackers_l1[cpu][index].conf);

  // update CS only if confidence is zero
  if (trackers_l1[cpu][index].conf == 0)
    trackers_l1[cpu][index].last_stride = stride;

  last_signature = trackers_l1[cpu][index].signature;
  // update complex stride(CPLX) confidence
  DPT_l1[cpu][last_signature].conf = update_conf(stride, DPT_l1[cpu][last_signature].delta, DPT_l1[cpu][last_signature].conf);

  // update CPLX only if confidence is zero
  if (DPT_l1[cpu][last_signature].conf == 0)
    DPT_l1[cpu][last_signature].delta = stride;

  // calculate and update new signature in IP table
  signature = update_sig_l1(last_signature, stride);
  trackers_l1[cpu][index].signature = signature;

  // check GHB for stream IP
  check_for_stream_l1(index, cl_addr, cpu);

  SIG_DP(cout << ip << ", " << cache_hit << ", " << cl_addr << ", " << addr << ", " << stride << "; ";
         cout << last_signature << ", " << DPT_l1[cpu][last_signature].delta << ", " << DPT_l1[cpu][last_signature].conf << "; ";
         cout << trackers_l1[cpu][index].last_stride << ", " << stride << ", " << trackers_l1[cpu][index].conf << ", "
              << "; ";);

  if (trackers_l1[cpu][index].str_valid == 1) { // stream IP
                                                // for stream, prefetch with twice the usual degree
    prefetch_degree = prefetch_degree * 2;
    for (int i = 0; i < prefetch_degree; i++) {
      uint64_t pf_address = 0;

      if (trackers_l1[cpu][index].str_dir == 1) { // +ve stream
        pf_address = (cl_addr + i + 1) << LOG2_BLOCK_SIZE;
        metadata = encode_metadata(1, S_TYPE, spec_nl[cpu]); // stride is 1
      } else {                                               // -ve stream
        pf_address = (cl_addr - i - 1) << LOG2_BLOCK_SIZE;
        metadata = encode_metadata(-1, S_TYPE, spec_nl[cpu]); // stride is -1
      }

      // Check if prefetch address is in same 4 KB page
      if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) {
        break;
      }

      prefetch_line(ip, addr, pf_address, FILL_L1, metadata);
      num_prefs++;
      SIG_DP(cout << "1, ");
    }

  } else if (trackers_l1[cpu][index].conf > 1 && trackers_l1[cpu][index].last_stride != 0) { // CS IP
    for (int i = 0; i < prefetch_degree; i++) {
      uint64_t pf_address = (cl_addr + (trackers_l1[cpu][index].last_stride * (i + 1))) << LOG2_BLOCK_SIZE;

      // Check if prefetch address is in same 4 KB page
      if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) {
        break;
      }

      metadata = encode_metadata(trackers_l1[cpu][index].last_stride, CS_TYPE, spec_nl[cpu]);
      prefetch_line(ip, addr, pf_address, FILL_L1, metadata);
      num_prefs++;
      SIG_DP(cout << trackers_l1[cpu][index].last_stride << ", ");
    }
  } else if (DPT_l1[cpu][signature].conf >= 0 && DPT_l1[cpu][signature].delta != 0) { // if conf>=0, continue looking for delta
    int pref_offset = 0, i = 0;                                                       // CPLX IP
    for (i = 0; i < prefetch_degree; i++) {
      pref_offset += DPT_l1[cpu][signature].delta;
      uint64_t pf_address = ((cl_addr + pref_offset) << LOG2_BLOCK_SIZE);

      // Check if prefetch address is in same 4 KB page
      if (((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) || (DPT_l1[cpu][signature].conf == -1) || (DPT_l1[cpu][signature].delta == 0)) {
        // if new entry in DPT or delta is zero, break
        break;
      }

      // we are not prefetching at L2 for CPLX type, so encode delta as 0
      metadata = encode_metadata(0, CPLX_TYPE, spec_nl[cpu]);
      if (DPT_l1[cpu][signature].conf > 0) { // prefetch only when conf>0 for CPLX
        prefetch_line(ip, addr, pf_address, FILL_L1, metadata);
        num_prefs++;
        SIG_DP(cout << pref_offset << ", ");
      }
      signature = update_sig_l1(signature, DPT_l1[cpu][signature].delta);
    }
  }

  // if no prefetches are issued till now, speculatively issue a next_line prefetch
  if (num_prefs == 0 && spec_nl[cpu] == 1) { // NL IP
    uint64_t pf_address = ((addr >> LOG2_BLOCK_SIZE) + 1) << LOG2_BLOCK_SIZE;
    metadata = encode_metadata(1, NL_TYPE, spec_nl[cpu]);
    prefetch_line(ip, addr, pf_address, FILL_L1, metadata);
    SIG_DP(cout << "1, ");
  }

  SIG_DP(cout << endl);

  // update the IP table entries
  trackers_l1[cpu][index].last_cl_offset = cl_offset;
  trackers_l1[cpu][index].last_page = curr_page;

  // update GHB
  // search for matching cl addr
  int ghb_index = 0;
  for (ghb_index = 0; ghb_index < NUM_GHB_ENTRIES; ghb_index++)
    if (cl_addr == ghb_l1[cpu][ghb_index])
      break;
  // only update the GHB upon finding a new cl address
  if (ghb_index == NUM_GHB_ENTRIES) {
    for (ghb_index = NUM_GHB_ENTRIES - 1; ghb_index > 0; ghb_index--)
      ghb_l1[cpu][ghb_index] = ghb_l1[cpu][ghb_index - 1];
    ghb_l1[cpu][0] = cl_addr;
  }

  return metadata_in;
}

uint32_t CACHE::prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{
  return metadata_in;
}
void CACHE::prefetcher_cycle_operate() {}

void CACHE::prefetcher_final_stats() {}
