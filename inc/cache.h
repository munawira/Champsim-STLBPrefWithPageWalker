#ifndef CACHE_H
#define CACHE_H

#include <functional>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "champsim.h"
#include "delay_queue.hpp"
#include "memory_class.h"
#include "ooo_cpu.h"
#include "operable.h"

// AADDED BY MUNAWIRA FOR MORRIGAN
//  PAGE
extern uint32_t PAGE_TABLE_LATENCY, SWAP_LATENCY;

#define P2TLB 0

// Free Prefetching

// flag for demand page walks
#define ENABLE_FP 1

// flag for prefetch page walks
#define ENABLE_PREF_FP 1

// Lookahead Depth (0:disable)
#define LA_DEPTH 0

// Replacement Policy for Markov's prediction table --> [ 0:LRU, 1:LFU, 2:RANDOM, 3:... ]
#define RP_MP 1

// Number of prediction table entries you randomly select from for eviction in LFU replacement policy for the Markov instruction TLB prefetcher
#define LLIMIT 5

// Number of bits for the confidence counters of the RP_SUC_MP replacement policy
#define CNF_BITS 3

// Number of successors of Markov I-TLB Prefetcher
#define SUCCESSORS 2

// Replacement Policy for the successors of Markov's prediction table --> [ 0:FIFO, 1:RANDOM, 2:CUSTOM, 3:... ]
#define RP_SUC_MP 2

// Number of STLB instruction misses to reset frequency field of Markov Prefetcher used for the LFU policy
#define RESET_FREQ 5000

#define PML4_SET 2
#define PDP_SET 4
#define PD_SET 16

#define PML4_WAY 1
#define PDP_WAY 1
#define PD_WAY 2

#define FCTB_SIZE 4

#define IS_STLB 2

// AADDED BY MUNAWIRA FOR MORRIGAN
//////////////////////////////////

// virtual address space prefetching
#define VA_PREFETCH_TRANSLATION_LATENCY 2

extern std::array<O3_CPU*, NUM_CPUS> ooo_cpu;

class CACHE : public champsim::operable, public MemoryRequestConsumer, public MemoryRequestProducer
{
public:
  uint32_t cpu;
  const std::string NAME;
  uint8_t cache_type;
  const uint32_t NUM_SET, NUM_WAY, WQ_SIZE, RQ_SIZE, PQ_SIZE, MSHR_SIZE;
  const uint32_t HIT_LATENCY, FILL_LATENCY, OFFSET_BITS;
  std::vector<BLOCK> block{NUM_SET * NUM_WAY};
  const uint32_t MAX_READ, MAX_WRITE;
  uint32_t reads_available_this_cycle, writes_available_this_cycle;
  const bool prefetch_as_load;
  const bool match_offset_bits;
  const bool virtual_prefetch;
  bool ever_seen_data = false;
  const unsigned pref_activate_mask = (1 << static_cast<int>(LOAD)) | (1 << static_cast<int>(PREFETCH));

  // prefetch stats
  uint64_t pf_requested = 0, pf_issued = 0, pf_useful = 0, pf_useless = 0, pf_fill = 0;

  // queues
  champsim::delay_queue<PACKET> RQ{RQ_SIZE, HIT_LATENCY}, // read queue
      PQ{PQ_SIZE, HIT_LATENCY},                           // prefetch queue
      VAPQ{PQ_SIZE, VA_PREFETCH_TRANSLATION_LATENCY},     // virtual address prefetch queue
      WQ{WQ_SIZE, HIT_LATENCY},                           // write queue
      PB{PQ_SIZE,HIT_LATENCY};

  std::list<PACKET> MSHR; // MSHR

  uint64_t sim_access[NUM_CPUS][NUM_TYPES] = {}, sim_hit[NUM_CPUS][NUM_TYPES] = {}, sim_miss[NUM_CPUS][NUM_TYPES] = {}, roi_access[NUM_CPUS][NUM_TYPES] = {},
           roi_hit[NUM_CPUS][NUM_TYPES] = {}, roi_miss[NUM_CPUS][NUM_TYPES] = {};

  uint64_t RQ_ACCESS = 0, RQ_MERGED = 0, RQ_FULL = 0, RQ_TO_CACHE = 0, PQ_ACCESS = 0, PQ_MERGED = 0, PQ_FULL = 0, PQ_TO_CACHE = 0, WQ_ACCESS = 0, WQ_MERGED = 0,
           WQ_FULL = 0, WQ_FORWARD = 0, WQ_TO_CACHE = 0;

  uint64_t total_miss_latency = 0;

  ////////////////////////////////////////////////
  // ADDED BY MUNAWIRA FOR MORRIGAN PREFETCHER
  uint64_t fctb[FCTB_SIZE][4];

  uint64_t timer;

  int selector;

  map<uint64_t, uint64_t> code_footprint;
  map<uint64_t, uint64_t>::iterator it;

  // prefetch stats
  uint64_t pf_hits_pq, pf_misses_pq,pf_hits_pb,pf_misses_pb, pf_swap, pf_dupli, pf_free, pf_real, previous_iva, previous_ip, fctb_hits, fctb_misses, hit_prefetches_lad,
      issued_prefetches_lad, pf_total_pq, morrigan_filter_hits, irip_hits, sdp_hits,
      bpbp[5],    // 0: # instruction prefetches 1: portion of instruction prefetches that are in the same page 2: portion of beyond page boundaries instruction
                  // prefetches 3: beyond page boundaries prefetches that hit in the TLB hierarchy 4: beyond page boundaries prefetches that miss in the TLB
      instr_miss, // Added by Munawira
      data_miss,  // Added by Munawira
      instr_trans_miss, // Added by Munawira
      data_trans_miss;  // Added by Munawira

  uint64_t free_distance_table[14];

  uint64_t pml4[PML4_SET][PML4_WAY], pdp[PDP_SET][PDP_WAY], pd[PD_SET][PD_WAY];
  uint64_t pml4_lru[PML4_SET][PML4_WAY], pdp_lru[PDP_SET][PDP_WAY], pd_lru[PD_SET][PD_WAY];
  uint64_t mmu_cache_demand_hits[4], mmu_cache_prefetch_hits[4];
  uint64_t mmu_timer;

  uint64_t pagetable_mr_hit_ratio[4][4];
  uint64_t rfhits[2];
  uint64_t free_hits[14];

  uint64_t decay_timer;
  uint64_t decay_conf_timer;

  uint64_t stlb_misses[2]; // 0 is data, 1 is instruction

  // functions
  pair<int, int> check_hit_stlb_pq(uint64_t vpn);

  // ADDED BY MUNAWIRA FOR MORRIGAN PREFETCHER
  ////////////////////////////////////////////////

  // functions
  int add_rq(PACKET* packet) override;
  int add_wq(PACKET* packet) override;
  int add_pq(PACKET* packet) override;

  void return_data(PACKET* packet) override;
  void operate() override;
  void operate_writes();
  void operate_reads();

  uint32_t get_occupancy(uint8_t queue_type, uint64_t address) override;
  uint32_t get_size(uint8_t queue_type, uint64_t address) override;

  uint32_t get_set(uint64_t address);
  uint32_t get_way(uint64_t address, uint32_t set);

  int invalidate_entry(uint64_t inval_addr);
  int prefetch_line(uint64_t pf_addr, bool fill_this_level, uint32_t prefetch_metadata);
  int prefetch_line(uint64_t ip, uint64_t base_addr, uint64_t pf_addr, bool fill_this_level, uint32_t prefetch_metadata); // deprecated

  void add_mshr(PACKET* packet);
  void va_translate_prefetches();

  void handle_fill();
  void handle_writeback();
  void handle_read();
  void handle_prefetch();

  void readlike_hit(std::size_t set, std::size_t way, PACKET& handle_pkt);
  bool readlike_miss(PACKET& handle_pkt);
  bool filllike_miss(std::size_t set, std::size_t way, PACKET& handle_pkt);

  bool should_activate_prefetcher(int type);

  void print_deadlock() override;

  ////////////////////////////////////////////////
  // ADDED BY MUNAWIRA FOR MORRIGAN PREFETCHER
  void stlb_prefetcher_initialize(), stlb_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr),
      stlb_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, int answer, int warmup, int*free_indexes, uint64_t instr_id,
                              int iflag),
      stlb_prefetcher_final_stats(uint64_t prefetches, uint64_t hits, uint64_t misses, uint64_t swap, uint64_t dupli, uint64_t free, uint64_t real,
                                  uint64_t*mmu_cache_demand_hits, uint64_t*mmu_cache_prefetch_hits, uint64_t*rfhits, uint64_t*free_hits, uint64_t mr[4][4],
                                  uint64_t stlb_misses[2]);

  int search_pml4(uint64_t address), search_pdp(uint64_t address), search_pd(uint64_t address);
  void lru_pml4(uint64_t timer, uint64_t address), lru_pdp(uint64_t timer, uint64_t address), lru_pd(uint64_t timer, uint64_t address);
  void print_pml4(), print_pdp(), print_pd();
  void free_prefetching(uint64_t ip, uint64_t addr, int cache_line_position_n, uint64_t pf_addr, int* free_indexes, uint64_t instr_id, int type, int iflag);

  // lookahead prefetching for Markov I-TLB Prefetcher
  void lookahead_prefetching(uint64_t ip, uint64_t addr, uint64_t pf_addr, uint64_t instr_id, int type, int iflag, int depth);

  // FCTB functions
  int fctb_replacement_policy(), search_fctb(uint64_t current_vpn);
  void print_fctb(), refresh_fctb(uint64_t current_cycle);

  int prefetch_page(uint64_t ip, uint64_t base_addr, uint64_t pf_addr, int fill_level, int pq_id, int free, int update_free, int free_distance, uint64_t id,
                    int type, int iflag, int lad, int confidence, int irip);
  int* sorted_free_distances();

  // ADDED BY MUNAWIRA FOR MORRIGAN PREFETCHER
  ////////////////////////////////////////////////

#include "cache_modules.inc"

  const repl_t repl_type;
  const pref_t pref_type;

  // constructor
  CACHE(std::string v1, double freq_scale, unsigned fill_level, uint32_t v2, int v3, uint32_t v5, uint32_t v6, uint32_t v7, uint32_t v8, uint32_t hit_lat,
        uint32_t fill_lat, uint32_t max_read, uint32_t max_write, std::size_t offset_bits, bool pref_load, bool wq_full_addr, bool va_pref,
        unsigned pref_act_mask, MemoryRequestConsumer* ll, pref_t pref, repl_t repl)
      : champsim::operable(freq_scale), MemoryRequestConsumer(fill_level), MemoryRequestProducer(ll), NAME(v1), NUM_SET(v2), NUM_WAY(v3), WQ_SIZE(v5),
        RQ_SIZE(v6), PQ_SIZE(v7), MSHR_SIZE(v8), HIT_LATENCY(hit_lat), FILL_LATENCY(fill_lat), OFFSET_BITS(offset_bits), MAX_READ(max_read),
        MAX_WRITE(max_write), prefetch_as_load(pref_load), match_offset_bits(wq_full_addr), virtual_prefetch(va_pref), pref_activate_mask(pref_act_mask),
        repl_type(repl), pref_type(pref)
  {

    ////////////////////////////////////////////////
    // ADDED BY MUNAWIRA FOR MORRIGAN PREFETCHER
    //  fctb initialization loop
    for (uint32_t j = 0; j < FCTB_SIZE; j++) {
      fctb[j][0] = 0;
      fctb[j][1] = 0;
      fctb[j][2] = 0;
      fctb[j][3] = 0;
    }

    for (int i = 0; i < PML4_SET; i++) {
      for (int j = 0; j < PML4_WAY; j++) {
        pml4[i][j] = 0;
        pml4_lru[i][j] = 0;
      }
    }

    for (int i = 0; i < PDP_SET; i++) {
      for (int j = 0; j < PDP_WAY; j++) {
        pdp[i][j] = 0;
        pdp_lru[i][j] = 0;
      }
    }

    for (int i = 0; i < PD_SET; i++) {
      for (int j = 0; j < PD_WAY; j++) {
        pd[i][j] = 0;
        pd_lru[i][j] = 0;
      }
    }

    mmu_timer = 0;

    for (int i = 0; i < 4; i++) {
      mmu_cache_demand_hits[i] = 0;
      mmu_cache_prefetch_hits[i] = 0;
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        pagetable_mr_hit_ratio[i][j] = 0;
      }
    }

    for (int i = 0; i < 14; i++) {
      free_hits[i] = 0;
    }

    rfhits[0] = 0;
    rfhits[1] = 0;
    decay_timer = 0;
    timer = 0;

    pf_hits_pq = 0;
    pf_misses_pq = 0;
    pf_swap = 0;
    pf_dupli = 0;
    pf_free = 0;
    pf_real = 0;
    previous_iva = 0;
    previous_ip = 0;
    stlb_misses[0] = 0;
    stlb_misses[1] = 0;
    fctb_hits = 0;
    fctb_misses = 0;
    issued_prefetches_lad = 0;
    hit_prefetches_lad = 0;
    pf_total_pq = 0;

    selector = 0;

    irip_hits = 0;
    sdp_hits = 0;

    bpbp[0] = 0;
    bpbp[1] = 0;
    bpbp[2] = 0;
    bpbp[3] = 0;
    bpbp[4] = 0;

    morrigan_filter_hits = 0;

    pf_requested = 0;
    pf_issued = 0;
    pf_useful = 0;
    pf_useless = 0;
    pf_fill = 0;

    cout << "CACHE NAME" << NAME<< endl;
    if (NAME == "cpu0_STLB")
      cache_type = IS_STLB;

    // ADDED BY MUNAWIRA FOR MORRIGAN PREFETCHER
    ////////////////////////////////////////////////
  }
};

#endif
