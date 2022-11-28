#include "cache.h"

#include <algorithm>
#include <iterator>

#include "champsim.h"
#include "champsim_constants.h"
#include "util.h"
#include "vmem.h"
#include "bits/stdc++.h"

#ifndef SANITY_CHECK
#define NDEBUG
#endif

extern VirtualMemory vmem;
extern uint8_t warmup_complete[NUM_CPUS];
int first_access =1;

////////////////////////////////////////////////
// ADDED below BY MUNAWIRA FOR MORRIGAN PREFETCHER

void CACHE::print_fctb()
{
  for (int i = 0; i < FCTB_SIZE; i++)
    cout << hex << fctb[i][0] << ", " << dec << fctb[i][1] << ", " << fctb[i][2] << ", " << fctb[i][3] << endl;
}

int CACHE::search_fctb(uint64_t current_vpn)
{
  int fctb_found_pos = -10, acc;
  for (int f = 0; f < FCTB_SIZE; f++) {
    if (fctb[f][0] == current_vpn)
      return f;

    for (int t = 1; t <= fctb[f][1]; t++) {
      if ((fctb[f][0] - t) == current_vpn)
        return f;
    }

    acc = 1;
    for (int t = fctb[f][1] + 1; t < 8; t++) {
      if ((fctb[f][0] + acc) == current_vpn)
        return f;
      acc++;
    }
  }
  return -10;
}

int CACHE::fctb_replacement_policy()
{
  uint64_t lru_min = fctb[0][2];
  int lru_victim = 0;
  for (int f = 1; f < FCTB_SIZE; f++) {
    if (fctb[f][2] < lru_min) {
      lru_min = fctb[f][2];
      lru_victim = f;
    }
  }
  return lru_victim;
}

void CACHE::refresh_fctb(uint64_t current_cycle)
{
  for (int f = 0; f < FCTB_SIZE; f++) {
    // if((current_cycle - fctb[f][2]) > PAGE_TABLE_LATENCY)
    if (current_cycle > 1 * (fctb[f][2] + fctb[f][3])) {
      fctb[f][0] = 0;
      fctb[f][1] = 0;
      fctb[f][2] = 0;
      fctb[f][3] = 0;
    }
  }
}

int* CACHE::sorted_free_distances()
{
  /*
     uint64_t sorted_table[14];

     for(int i=0; i<14; i++)
     sorted_table[i] = free_distance_table[i];

     int n = sizeof(sorted_table)/sizeof(sorted_table[0]);

     sort(sorted_table, sorted_table+n);

  //for(int i=0; i<14; i++){
  //    cout << sorted_table[i] << ", ";
  //}
  //cout << endl;
  */

  /* pair sort */
  int n = sizeof(free_distance_table) / sizeof(free_distance_table[0]);
  pair<uint64_t, int>* pairf;
  pairf = (pair<uint64_t, int>*)malloc(14 * sizeof(pair<uint64_t, int>));

  for (int i = 0; i < 14; i++) {
    pairf[i].second = i - 6;
    if (i <= 6)
      pairf[i].second--;
    pairf[i].first = free_distance_table[i];
  }

  sort(pairf, pairf + n);

  /*
     for(int i=0; i<14; i++){
     cout << pairf[i].first << ", ";
     }
     cout << endl;
     for(int i=0; i<14; i++){
     cout << pairf[i].second << ", ";
     }
     */

  int* indexes;
  indexes = (int*)malloc(14 * sizeof(int));
  for (int i = 0; i < 14; i++) {
    indexes[i] = pairf[i].second;
  }

  free(pairf);

  return indexes;

  // return pairf;
}

void CACHE::print_pml4()
{
  for (int s = 0; s < PML4_SET; s++) {
    for (int w = 0; w < PML4_WAY; w++)
      cout << pml4[s][w] << ", " << pml4_lru[s][w];
    cout << endl;
  }
}

void CACHE::print_pdp()
{
  for (int s = 0; s < PDP_SET; s++) {
    for (int w = 0; w < PDP_WAY; w++)
      cout << pdp[s][w] << ", " << pdp_lru[s][w];
    cout << endl;
  }
}

void CACHE::print_pd()
{
  for (int s = 0; s < PD_SET; s++) {
    for (int w = 0; w < PD_WAY; w++)
      cout << pd[s][w] << ", " << pd_lru[s][w] << " ;; ";
    cout << endl;
  }
}

int CACHE::search_pml4(uint64_t address)
{
  for (int s = 0; s < PML4_SET; s++) {
    for (int w = 0; w < PML4_WAY; w++) {
      if (pml4[s][w] == address)
        return 1;
    }
  }
  return 0;
}

int CACHE::search_pdp(uint64_t address)
{
  for (int s = 0; s < PDP_SET; s++) {
    for (int w = 0; w < PDP_WAY; w++) {
      if (pdp[s][w] == address)
        return 1;
    }
  }
  return 0;
}

int CACHE::search_pd(uint64_t address)
{
  for (int s = 0; s < PD_SET; s++) {
    for (int w = 0; w < PD_WAY; w++) {
      if (pd[s][w] == address)
        return 1;
    }
  }
  return 0;
}

void CACHE::lru_pml4(uint64_t timer, uint64_t address)
{
  for (int s = 0; s < PML4_SET; s++) {
    for (int w = 0; w < PML4_WAY; w++) {
      if (pml4[s][w] == address) {
        pml4_lru[s][w] = timer;
        return;
      }
    }
  }

  uint64_t min = pml4_lru[0][0];
  int victim[2] = {0, 0};

  for (int s = 0; s < PML4_SET; s++) {
    for (int w = 0; w < PML4_WAY; w++) {
      if (pml4_lru[s][w] < min) {
        min = pml4_lru[s][w];
        victim[0] = s;
        victim[1] = w;
      }
    }
  }

  pml4[victim[0]][victim[1]] = address;
  pml4_lru[victim[0]][victim[1]] = timer;
  return;
}

void CACHE::lru_pdp(uint64_t timer, uint64_t address)
{
  for (int s = 0; s < PDP_SET; s++) {
    for (int w = 0; w < PDP_WAY; w++) {
      if (pdp[s][w] == address) {
        pdp_lru[s][w] = timer;
        return;
      }
    }
  }

  uint64_t min = pdp_lru[0][0];
  int victim[2] = {0, 0};

  for (int s = 0; s < PDP_SET; s++) {
    for (int w = 0; w < PDP_WAY; w++) {
      if (pdp_lru[s][w] < min) {
        min = pdp_lru[s][w];
        victim[0] = s;
        victim[1] = w;
      }
    }
  }

  pdp[victim[0]][victim[1]] = address;
  pdp_lru[victim[0]][victim[1]] = timer;
  return;
}
void CACHE::lru_pd(uint64_t timer, uint64_t address)
{
  int bits = ceil(log2(PD_SET));
  int index = address & ((1 << bits) - 1);
  // cout << "Index: " << index << endl;

  for (int w = 0; w < PD_WAY; w++) {
    if (pd[index][w] == address) {
      pd_lru[index][w] = timer;
      return;
    }
  }

  uint64_t min = pd_lru[index][0];
  int victim[2] = {index, 0};

  for (int w = 0; w < PD_WAY; w++) {
    if (pd_lru[index][w] < min) {
      min = pd_lru[index][w];
      victim[0] = index;
      victim[1] = w;
    }
  }

  pd[victim[0]][victim[1]] = address;
  pd_lru[victim[0]][victim[1]] = timer;

  return;
}

int compare(const void* a, const void* b) { return (*(int*)a - *(int*)b); }

uint64_t l2pf_access = 0;


int CACHE::prefetch_page(uint64_t ip, uint64_t base_addr, uint64_t pf_addr, int fill_level, int pq_id, int free, int update_free, int free_distance,
                         uint64_t id, int type, int iflag, int lad, int confidence, int irip)
{
  int index = 0, debug = 0, flag = 0, fctb_search = -10;
  bool vapq_full = false;
  //uint64_t temp = va_to_pa_prefetch(cpu, base_addr, pf_addr), foo;
  

  if(!free)
  	fctb_search = search_fctb(pf_addr);

  if(pq_id == 0){
  	pf_requested++;
  	pf_total_pq++;
  }

  // if(pq_id == 0){
  // 	auto vapq_entry = std::find_if(PQ.begin(), VAPQ.end(), eq_addr<PACKET>(pf_addr, OFFSET_BITS));
  //   vapq_full = (VAPQ.full());
  //   if (vapq_entry == VAPQ.end())
  //     index = -1;
  // }
  // else{
  // 	cout << "I am using only one PQ" << endl;
  // }

  // if(vapq_full){
  //   //cout << "VAPQ FULL : " << VAPQ.size() << ": " << PQ_SIZE << endl;
  //   return 0;
  // }

  // //cout << "VAPQ NOT FULL" << endl;
  // if(pq_id != 2){
  //   if(index != -1){
  //    // cout << "Duplicate" << endl;
  //     if(debug)
  //       cout << "Duplicate in the Prefetch Queue: " << pf_addr << endl;
  //     return 0;
  //   }
  // }

   	PACKET pf_packet;


  	pf_packet.fill_level = fill_level;
  	pf_packet.cpu = cpu;
  	pf_packet.address = pf_addr;
  	pf_packet.v_address = pf_addr;
  	pf_packet.ip = ip;
  	pf_packet.type = TRANSLATION;
  	pf_packet.event_cycle = current_core_cycle[cpu];
  	pf_packet.free_bit = free;
  	pf_packet.free_distance = free_distance;
  	pf_packet.lad = lad;
  	pf_packet.conf = confidence;
  	pf_packet.irip = irip;
    pf_packet.is_instr_addr =1;
    pf_packet.instruction =1;
    pf_packet.is_stlb_prefetch =1;


    //pf_packet.to_return = {this};

    //auto it = MSHR.insert(std::end(MSHR), pf_packet);
    // it->cycle_enqueued = current_cycle;
    // it->event_cycle = std::numeric_limits<uint64_t>::max();


  	if(fctb_search == -10){
  		fctb_misses++;
  		pf_packet.event_cycle = current_core_cycle[cpu];
  		pf_packet.free_bit = free;
  	}
  	else{
  		fctb_hits++;
  		pf_packet.event_cycle = fctb[fctb_search][2];
  		pf_packet.free_bit = 1;
  	}

  	if(free){
  		if(lad == 0)
  			pf_free++;
  	}
  	else{
  		if(fctb_search != -10){
  			pf_free++;
  		}
  		else{
  			pf_real++;

  			int victim_entry = fctb_replacement_policy();
  			fctb[victim_entry][0] = pf_addr;
  			fctb[victim_entry][1] = (pf_addr & 0x07);
  			fctb[victim_entry][2] = current_core_cycle[cpu];
  			
  		}
  	}

    if(STLB_PB.size() != STLB_PB_SIZE){
      

      auto it = STLB_PB.insert(std::end(STLB_PB), pf_packet);
      stlb_pb_added++;
      //cout << "STLB PB INSERTED " << pf_packet.address << endl; 
              
    }
 
  	if(pq_id == 0){

      //cout<< "ADD PQ " <<endl;
      lower_level->add_pq(&pf_packet);
     //add_pq(&pf_packet);
     // VAPQ.push_back(pf_packet);
    }else
  		cout << "I am using only one PQ" << endl;

  	if(lad == 1)
  		issued_prefetches_lad++;

  	pf_issued++;
  	return 1;
}

// ADDED BY MUNAWIRA FOR MORRIGAN PREFETCHER
////////////////////////////////////////////////

void CACHE::handle_fill()
{
  while (writes_available_this_cycle > 0) {
    auto fill_mshr = MSHR.begin();
    if (fill_mshr == std::end(MSHR) || fill_mshr->event_cycle > current_cycle)
      return;
    
    // if((fill_mshr->is_stlb_prefetch == 1) && (NAME == "cpu0_STLB")){
    //   cout << "FILLING FOR STLB PREFETCH" << endl;
      
    //   auto it = STLB_PB.insert(std::end(STLB_PB),*fill_mshr);
    //   stlb_pb_added++;

    // }else{  

    // find victim
    uint32_t set = get_set(fill_mshr->address);

    auto set_begin = std::next(std::begin(block), set * NUM_WAY);
    auto set_end = std::next(set_begin, NUM_WAY);
    auto first_inv = std::find_if_not(set_begin, set_end, is_valid<BLOCK>());
    uint32_t way = std::distance(set_begin, first_inv);
    if (way == NUM_WAY)
      way = impl_replacement_find_victim(fill_mshr->cpu, fill_mshr->instr_id, set, &block.data()[set * NUM_WAY], fill_mshr->ip, fill_mshr->address,
                                         fill_mshr->type);

    bool success = filllike_miss(set, way, *fill_mshr);
    if (!success)
      return;

    if (way != NUM_WAY) {
      // update processed packets
      fill_mshr->data = block[set * NUM_WAY + way].data;

      for (auto ret : fill_mshr->to_return)
        ret->return_data(&(*fill_mshr));
    }
    
    MSHR.erase(fill_mshr);
    writes_available_this_cycle--;
  }
}

void CACHE::handle_writeback()
{
  while (writes_available_this_cycle > 0) {
    if (!WQ.has_ready())
      return;

    // handle the oldest entry
    PACKET& handle_pkt = WQ.front();

    // access cache
    uint32_t set = get_set(handle_pkt.address);
    uint32_t way = get_way(handle_pkt.address, set);

    BLOCK& fill_block = block[set * NUM_WAY + way];

    if (way < NUM_WAY) // HIT
    {
      impl_replacement_update_state(handle_pkt.cpu, set, way, fill_block.address, handle_pkt.ip, 0, handle_pkt.type, 1);

      // COLLECT STATS
      sim_hit[handle_pkt.cpu][handle_pkt.type]++;
      sim_access[handle_pkt.cpu][handle_pkt.type]++;

      // mark dirty
      fill_block.dirty = 1;
    } else // MISS
    {
      bool success;
      if (handle_pkt.type == RFO && handle_pkt.to_return.empty()) {
        success = readlike_miss(handle_pkt);
      } else {
        // find victim
        auto set_begin = std::next(std::begin(block), set * NUM_WAY);
        auto set_end = std::next(set_begin, NUM_WAY);
        auto first_inv = std::find_if_not(set_begin, set_end, is_valid<BLOCK>());
        way = std::distance(set_begin, first_inv);
        if (way == NUM_WAY)
          way = impl_replacement_find_victim(handle_pkt.cpu, handle_pkt.instr_id, set, &block.data()[set * NUM_WAY], handle_pkt.ip, handle_pkt.address,
                                             handle_pkt.type);

        success = filllike_miss(set, way, handle_pkt);
      }

      if (!success)
        return;
    }

    // remove this entry from WQ
    writes_available_this_cycle--;
    WQ.pop_front();
  }
}

void CACHE::handle_read()
{
  while (reads_available_this_cycle > 0) {

    if (!RQ.has_ready())
      return;

    // handle the oldest entry
    PACKET& handle_pkt = RQ.front();

    // A (hopefully temporary) hack to know whether to send the evicted paddr or
    // vaddr to the prefetcher
    ever_seen_data |= (handle_pkt.v_address != handle_pkt.ip);

    uint32_t set = get_set(handle_pkt.address);
    uint32_t way = get_way(handle_pkt.address, set);

    if (way < NUM_WAY) // HIT
    {
      readlike_hit(set, way, handle_pkt);
    } else {
      bool success = readlike_miss(handle_pkt);
      if (!success)
        return;
    }

    // remove this entry from RQ
    RQ.pop_front();
    reads_available_this_cycle--;
  }
}

void CACHE::handle_prefetch()
{
  while (reads_available_this_cycle > 0) {
    if (!PQ.has_ready())
      return;
    
    // handle the oldest entry
    PACKET& handle_pkt = PQ.front();

    uint32_t set = get_set(handle_pkt.address);
    uint32_t way = get_way(handle_pkt.address, set);

    if (way < NUM_WAY) // HIT
    {
      if(handle_pkt.is_stlb_prefetch)
          pf_hits_pq++;
      readlike_hit(set, way, handle_pkt);
    } else {
      if(handle_pkt.is_stlb_prefetch)
          pf_misses_pq++;
      bool success = readlike_miss(handle_pkt);
      if (!success)
        return;
    }

    // remove this entry from PQ
    PQ.pop_front();
    reads_available_this_cycle--;
  }
}

void CACHE::readlike_hit(std::size_t set, std::size_t way, PACKET& handle_pkt)
{
  DP(if (warmup_complete[handle_pkt.cpu]) {
    std::cout << "[" << NAME << "] " << __func__ << " hit";
    std::cout << " instr_id: " << handle_pkt.instr_id << " address: " << std::hex << (handle_pkt.address >> OFFSET_BITS);
    std::cout << " full_addr: " << handle_pkt.address;
    std::cout << " full_v_addr: " << handle_pkt.v_address << std::dec;
    std::cout << " type: " << +handle_pkt.type;
    std::cout << " cycle: " << current_cycle << std::endl;
  });

  BLOCK& hit_block = block[set * NUM_WAY + way];

  handle_pkt.data = hit_block.data;

  // update prefetcher on load instruction //MUNA QUESTION: ANYTHING TO DO HERE
  if (should_activate_prefetcher(handle_pkt.type) && handle_pkt.pf_origin_level < fill_level) {
    cpu = handle_pkt.cpu;
    uint64_t pf_base_addr = (virtual_prefetch ? handle_pkt.v_address : handle_pkt.address) & ~bitmask(match_offset_bits ? 0 : OFFSET_BITS);
    handle_pkt.pf_metadata = impl_prefetcher_cache_operate(pf_base_addr, handle_pkt.ip, 1, handle_pkt.type, handle_pkt.pf_metadata);
  }

  // update replacement policy
  impl_replacement_update_state(handle_pkt.cpu, set, way, hit_block.address, handle_pkt.ip, 0, handle_pkt.type, 1);

  // COLLECT STATS
  sim_hit[handle_pkt.cpu][handle_pkt.type]++;
  sim_access[handle_pkt.cpu][handle_pkt.type]++;

  for (auto ret : handle_pkt.to_return)
    ret->return_data(&handle_pkt); // MUNA: Return data from Cache

  // update prefetch stats and reset prefetch bit
  if (hit_block.prefetch) {
    pf_useful++;
    hit_block.prefetch = 0;
  }
}

bool CACHE::filllike_miss(std::size_t set, std::size_t way, PACKET& handle_pkt)
{
  DP(if (warmup_complete[handle_pkt.cpu]) {
    std::cout << "[" << NAME << "] " << __func__ << " miss";
    std::cout << " instr_id: " << handle_pkt.instr_id << " address: " << std::hex << (handle_pkt.address >> OFFSET_BITS);
    std::cout << " full_addr: " << handle_pkt.address;
    std::cout << " full_v_addr: " << handle_pkt.v_address << std::dec;
    std::cout << " type: " << +handle_pkt.type;
    std::cout << " cycle: " << current_cycle << std::endl;
  });

  bool bypass = (way == NUM_WAY);
#ifndef LLC_BYPASS
  assert(!bypass);
#endif
  assert(handle_pkt.type != WRITEBACK || !bypass);

  BLOCK& fill_block = block[set * NUM_WAY + way];
  bool evicting_dirty = !bypass && (lower_level != NULL) && fill_block.dirty;
  uint64_t evicting_address = 0;

  if (!bypass) {
    if (evicting_dirty) {
      PACKET writeback_packet;

      writeback_packet.fill_level = lower_level->fill_level;
      writeback_packet.cpu = handle_pkt.cpu;
      writeback_packet.address = fill_block.address;
      writeback_packet.data = fill_block.data;
      writeback_packet.instr_id = handle_pkt.instr_id;
      writeback_packet.ip = 0;
      writeback_packet.type = WRITEBACK;

      auto result = lower_level->add_wq(&writeback_packet);
      if (result == -2)
        return false;
    }

    if (ever_seen_data)
      evicting_address = fill_block.address & ~bitmask(match_offset_bits ? 0 : OFFSET_BITS);
    else
      evicting_address = fill_block.v_address & ~bitmask(match_offset_bits ? 0 : OFFSET_BITS);

    if (fill_block.prefetch)
      pf_useless++;

    if (handle_pkt.type == PREFETCH)
      pf_fill++;

    fill_block.valid = true;
    fill_block.prefetch = (handle_pkt.type == PREFETCH && handle_pkt.pf_origin_level == fill_level);
    fill_block.dirty = (handle_pkt.type == WRITEBACK || (handle_pkt.type == RFO && handle_pkt.to_return.empty()));
    fill_block.address = handle_pkt.address;
    fill_block.v_address = handle_pkt.v_address;
    fill_block.data = handle_pkt.data;
    fill_block.ip = handle_pkt.ip;
    fill_block.cpu = handle_pkt.cpu;
    fill_block.instr_id = handle_pkt.instr_id;
  }

  if (warmup_complete[handle_pkt.cpu] && (handle_pkt.cycle_enqueued != 0))
    total_miss_latency += current_cycle - handle_pkt.cycle_enqueued;

  // update prefetcher
  cpu = handle_pkt.cpu;
  handle_pkt.pf_metadata =
      impl_prefetcher_cache_fill((virtual_prefetch ? handle_pkt.v_address : handle_pkt.address) & ~bitmask(match_offset_bits ? 0 : OFFSET_BITS), set, way,
                                 handle_pkt.type == PREFETCH, evicting_address, handle_pkt.pf_metadata);

  // update replacement policy
  impl_replacement_update_state(handle_pkt.cpu, set, way, handle_pkt.address, handle_pkt.ip, 0, handle_pkt.type, 0);

  // COLLECT STATS
  sim_miss[handle_pkt.cpu][handle_pkt.type]++;
  sim_access[handle_pkt.cpu][handle_pkt.type]++;

  return true;
}

bool CACHE::readlike_miss(PACKET& handle_pkt)
{
  DP(if (warmup_complete[handle_pkt.cpu]) {
    std::cout << "[" << NAME << "] " << __func__ << " miss";
    std::cout << " instr_id: " << handle_pkt.instr_id << " address: " << std::hex << (handle_pkt.address >> OFFSET_BITS);
    std::cout << " full_addr: " << handle_pkt.address;
    std::cout << " full_v_addr: " << handle_pkt.v_address << std::dec;
    std::cout << " type: " << +handle_pkt.type;
    std::cout << " cycle: " << current_cycle << std::endl;
  });

  bool morrigan_engage = false;
  // check mshr
  auto mshr_entry = std::find_if(MSHR.begin(), MSHR.end(), eq_addr<PACKET>(handle_pkt.address, OFFSET_BITS));
  bool mshr_full = (MSHR.size() == MSHR_SIZE);

  if (mshr_entry != MSHR.end()) // miss already inflight //MUNA TODO: Check if prefetch needed if miss already inflight
  {
    // update fill location
    mshr_entry->fill_level = std::min(mshr_entry->fill_level, handle_pkt.fill_level);

    packet_dep_merge(mshr_entry->lq_index_depend_on_me, handle_pkt.lq_index_depend_on_me);
    packet_dep_merge(mshr_entry->sq_index_depend_on_me, handle_pkt.sq_index_depend_on_me);
    packet_dep_merge(mshr_entry->instr_depend_on_me, handle_pkt.instr_depend_on_me);
    packet_dep_merge(mshr_entry->to_return, handle_pkt.to_return);

    if (mshr_entry->type == PREFETCH && handle_pkt.type != PREFETCH) {
      // Mark the prefetch as useful
      if (mshr_entry->pf_origin_level == fill_level)
        pf_useful++;

      uint64_t prior_event_cycle = mshr_entry->event_cycle;
      *mshr_entry = handle_pkt;

      // in case request is already returned, we should keep event_cycle
      mshr_entry->event_cycle = prior_event_cycle;
    }
  } else {
    if (mshr_full)  // not enough MSHR resource
      return false; // TODO should we allow prefetches anyway if they will not
                    // be filled to this level?

    // MUNA: BELOW: if entry not already found in MSHR and mshr is not full
    bool is_read = prefetch_as_load || (handle_pkt.type != PREFETCH); // MUNA: Check if morrigan prefetches as load

    // check to make sure the lower level queue has room for this read miss
    int queue_type = (is_read) ? 1 : 3;
    if (lower_level->get_occupancy(queue_type, handle_pkt.address) == lower_level->get_size(queue_type, handle_pkt.address))
      return false;

    // Allocate an MSHR
    if (handle_pkt.fill_level <= fill_level) {
      auto it = MSHR.insert(std::end(MSHR), handle_pkt);
      it->cycle_enqueued = current_cycle;
      it->event_cycle = std::numeric_limits<uint64_t>::max();
    }

    if (handle_pkt.fill_level <= fill_level) {
      handle_pkt.to_return = {this};
    } else {
      handle_pkt.to_return.clear();
    }

    if (!is_read) {
      lower_level->add_pq(&handle_pkt);
   
    } else {
 
      if ((this->NAME == "cpu0_STLB") && MORRIGAN && (handle_pkt.instruction!=0))
        morrigan_engage = true;


      if(!morrigan_engage){
         lower_level->add_rq(&handle_pkt);
      } else { // Only if Instruction STLB Miss and Morrigan Activated, This part of the code will engage
      //Handle to call only if DSTLB      // MUNA: BELOW ADDED BY MUNAWIRA FOR MORRIGAN PREFETCHER
     
        
        uint64_t pa, current_vpn = handle_pkt.v_address;
        int bits, rowhit = -1, victim, iflag = 0;
        pair<int, int> answer;

        int* free_indexes;
        free_indexes = sorted_free_distances();

        if (unsigned(handle_pkt.instruction) != 0) {
          stlb_misses[1]++;
          iflag = 1;
          decay_timer++;
          decay_conf_timer++;
          instr_trans_miss++;
          handle_pkt.is_instr_addr = 1;
          handle_pkt.is_data_addr = 0;
        } else {
          stlb_misses[0]++;
          data_trans_miss++;

          handle_pkt.is_instr_addr = 0;
          handle_pkt.is_data_addr = 1;
        }

        int fctb_found_pos = -10;

        if (iflag == 1) {
          refresh_fctb(current_cycle);
          fctb_found_pos = search_fctb(handle_pkt.v_address);

          //pa =vmem.va_to_pa(0,handle_pkt.v_address).first;
          //Search PQ replaced with this
        
          auto PB_entry = std::find_if(STLB_PB.begin(), STLB_PB.end(), eq_addr<PACKET>(handle_pkt.address, OFFSET_BITS));
          //bool PB_full = (STLB_PB.size() == STLB_PB_SIZE);
          answer = make_pair(1, 0);
          if (PB_entry == STLB_PB.end()){
            answer = make_pair(-1, -1);
          }      
          
        } else {
          answer = make_pair(-1, -1);
        }

        pair<uint64_t, uint64_t> v2p;
        if (answer.first == -1) {// MUNA: Not found in PB
        //cout << "NOT FOUND IN PB" << endl;
          if (iflag) {
            if (fctb_found_pos == -10) {
              int victim_entry = fctb_replacement_policy();
              fctb[victim_entry][0] = current_vpn;
              fctb[victim_entry][1] = (handle_pkt.v_address & 0x7000) / 4096;
              fctb[victim_entry][2] = current_cycle;
              fctb[victim_entry][3] = v2p.second;
            }
          }

          lower_level->add_rq(&handle_pkt);
          if (iflag == 1) {
            pf_misses_pb++;
          }


        } else { // Update stats in Prefetch queue, if data found in PQ
          //cout << "FOUND IN PB" << endl;
          pa =vmem.va_to_pa(0,handle_pkt.v_address).first;
          
        
          auto PB_entry = std::find_if(STLB_PB.begin(), STLB_PB.end(), eq_addr<PACKET>(handle_pkt.address, OFFSET_BITS));
          if (PB_entry->free_bit == 1 && PB_entry->free_distance != 0) {
            rfhits[1]++;
            free_hits[PB_entry->free_distance + 6 + (PB_entry->free_distance < 0) * 1]++;
          } else
            rfhits[0]++;

          if (warmup_complete[cpu]) {
            uint64_t num_cycles;

            if (PB_entry->irip) // MUNA: updated in STLB prefetcher operate
              irip_hits++;
            else
              sdp_hits++;
          }

          pf_hits_pb++;
          // cout << "Found in PB" << endl;

          handle_pkt.data = pa >> LOG2_PAGE_SIZE;
          handle_pkt.to_return = {this};
          return_data(&handle_pkt);
          
          STLB_PB.erase(PB_entry);



        }
        if (iflag == 1) {
          //cout<<"We're in STLB MISS 5" << endl;
          free_indexes = sorted_free_distances();
          stlb_prefetcher_operate(handle_pkt.v_address, handle_pkt.ip, 0, handle_pkt.type, answer.first, warmup_complete[cpu], free_indexes,
                                  handle_pkt.instr_id, iflag);
          // stlb_prefetcher_cache_fill(handle_pkt.v_address, 0, 0, 0, 0); //MUNA: No implementation specified
        }

        // MUNA: ABOVE ADDED BY MUNAWIRA  FOR MORRIGAN PREFETCHER
      }
    }
  }
  // update prefetcher on load instructions and prefetches from upper levels
  if (should_activate_prefetcher(handle_pkt.type) && handle_pkt.pf_origin_level < fill_level) {
    cpu = handle_pkt.cpu;
    uint64_t pf_base_addr = (virtual_prefetch ? handle_pkt.v_address : handle_pkt.address) & ~bitmask(match_offset_bits ? 0 : OFFSET_BITS);
    handle_pkt.pf_metadata = impl_prefetcher_cache_operate(pf_base_addr, handle_pkt.ip, 0, handle_pkt.type, handle_pkt.pf_metadata);
  }

  return true;
}

void CACHE::operate()
{
  operate_writes();
  operate_reads();

  impl_prefetcher_cycle_operate();
}

void CACHE::operate_writes()
{
  // perform all writes
  writes_available_this_cycle = MAX_WRITE;
  handle_fill();
  handle_writeback();

  WQ.operate();
}

void CACHE::operate_reads()
{
  // perform all reads
  reads_available_this_cycle = MAX_READ;
  handle_read();
  va_translate_prefetches();
  handle_prefetch();

  RQ.operate();
  PQ.operate();
  VAPQ.operate();
}

uint32_t CACHE::get_set(uint64_t address) { return ((address >> OFFSET_BITS) & bitmask(lg2(NUM_SET))); }

uint32_t CACHE::get_way(uint64_t address, uint32_t set)
{
  auto begin = std::next(block.begin(), set * NUM_WAY);
  auto end = std::next(begin, NUM_WAY);
  return std::distance(begin, std::find_if(begin, end, eq_addr<BLOCK>(address, OFFSET_BITS)));
}

int CACHE::invalidate_entry(uint64_t inval_addr)
{
  uint32_t set = get_set(inval_addr);
  uint32_t way = get_way(inval_addr, set);

  if (way < NUM_WAY)
    block[set * NUM_WAY + way].valid = 0;

  return way;
}

int CACHE::add_rq(PACKET* packet)
{
  assert(packet->address != 0);
  RQ_ACCESS++;

  DP(if (warmup_complete[packet->cpu]) {
    std::cout << "[" << NAME << "_RQ] " << __func__ << " instr_id: " << packet->instr_id << " address: " << std::hex << (packet->address >> OFFSET_BITS);
    std::cout << " full_addr: " << packet->address << " v_address: " << packet->v_address << std::dec << " type: " << +packet->type
              << " occupancy: " << RQ.occupancy();
  })

  // check for the latest writebacks in the write queue
  champsim::delay_queue<PACKET>::iterator found_wq = std::find_if(WQ.begin(), WQ.end(), eq_addr<PACKET>(packet->address, match_offset_bits ? 0 : OFFSET_BITS));

  if (found_wq != WQ.end()) {

    DP(if (warmup_complete[packet->cpu]) std::cout << " MERGED_WQ" << std::endl;)

    packet->data = found_wq->data;
    for (auto ret : packet->to_return)
      ret->return_data(packet);

    WQ_FORWARD++;
    return -1;
  }

  // check for duplicates in the read queue
  auto found_rq = std::find_if(RQ.begin(), RQ.end(), eq_addr<PACKET>(packet->address, OFFSET_BITS));
  if (found_rq != RQ.end()) {

    DP(if (warmup_complete[packet->cpu]) std::cout << " MERGED_RQ" << std::endl;)

    packet_dep_merge(found_rq->lq_index_depend_on_me, packet->lq_index_depend_on_me);
    packet_dep_merge(found_rq->sq_index_depend_on_me, packet->sq_index_depend_on_me);
    packet_dep_merge(found_rq->instr_depend_on_me, packet->instr_depend_on_me);
    packet_dep_merge(found_rq->to_return, packet->to_return);

    RQ_MERGED++;

    return 0; // merged index
  }

  // check occupancy
  if (RQ.full()) {
    RQ_FULL++;

    DP(if (warmup_complete[packet->cpu]) std::cout << " FULL" << std::endl;)

    return -2; // cannot handle this request
  }

  // if there is no duplicate, add it to RQ
  if (warmup_complete[cpu])
    RQ.push_back(*packet);
  else
    RQ.push_back_ready(*packet);

  DP(if (warmup_complete[packet->cpu]) std::cout << " ADDED" << std::endl;)

  RQ_TO_CACHE++;
  return RQ.occupancy();
}

int CACHE::add_wq(PACKET* packet)
{
  WQ_ACCESS++;

  DP(if (warmup_complete[packet->cpu]) {
    std::cout << "[" << NAME << "_WQ] " << __func__ << " instr_id: " << packet->instr_id << " address: " << std::hex << (packet->address >> OFFSET_BITS);
    std::cout << " full_addr: " << packet->address << " v_address: " << packet->v_address << std::dec << " type: " << +packet->type
              << " occupancy: " << RQ.occupancy();
  })

  // check for duplicates in the write queue
  champsim::delay_queue<PACKET>::iterator found_wq = std::find_if(WQ.begin(), WQ.end(), eq_addr<PACKET>(packet->address, match_offset_bits ? 0 : OFFSET_BITS));

  if (found_wq != WQ.end()) {

    DP(if (warmup_complete[packet->cpu]) std::cout << " MERGED" << std::endl;)

    WQ_MERGED++;
    return 0; // merged index
  }

  // Check for room in the queue
  if (WQ.full()) {
    DP(if (warmup_complete[packet->cpu]) std::cout << " FULL" << std::endl;)

    ++WQ_FULL;
    return -2;
  }

  // if there is no duplicate, add it to the write queue
  if (warmup_complete[cpu]) {
    WQ.push_back(*packet);
  } else
    WQ.push_back_ready(*packet);

  DP(if (warmup_complete[packet->cpu]) std::cout << " ADDED" << std::endl;)

  WQ_TO_CACHE++;
  WQ_ACCESS++;

  return WQ.occupancy();
}

int CACHE::prefetch_line(uint64_t pf_addr, bool fill_this_level, uint32_t prefetch_metadata)
{
  pf_requested++;

  PACKET pf_packet;
  pf_packet.type = PREFETCH;
  pf_packet.fill_level = (fill_this_level ? fill_level : lower_level->fill_level);
  pf_packet.pf_origin_level = fill_level;
  pf_packet.pf_metadata = prefetch_metadata;
  pf_packet.cpu = cpu;
  pf_packet.address = pf_addr;
  pf_packet.v_address = virtual_prefetch ? pf_addr : 0;

  if (virtual_prefetch) {
    if (!VAPQ.full()) {
      VAPQ.push_back(pf_packet);
      return 1;
    }
  } else {
    int result = add_pq(&pf_packet);
    if (result != -2) {
      if (result > 0) {
        pf_issued++;
      }
      return 1;
    }
  }

  return 0;
}

int CACHE::prefetch_line(uint64_t ip, uint64_t base_addr, uint64_t pf_addr, bool fill_this_level, uint32_t prefetch_metadata)
{
  static bool deprecate_printed = false;
  if (!deprecate_printed) {
    std::cout << "WARNING: The extended signature CACHE::prefetch_line(ip, "
                 "base_addr, pf_addr, fill_this_level, prefetch_metadata) is "
                 "deprecated."
              << std::endl;
    std::cout << "WARNING: Use CACHE::prefetch_line(pf_addr, fill_this_level, "
                 "prefetch_metadata) instead."
              << std::endl;
    deprecate_printed = true;
  }
  return prefetch_line(pf_addr, fill_this_level, prefetch_metadata);
}

void CACHE::va_translate_prefetches()
{
  int result = -2;
  if (this->NAME == "cpu0_STLB" && MORRIGAN) { // Page Walker goes through read queue
    if (!VAPQ.empty()){
      cout << "VAPQ Address:" << VAPQ.front().address;
      result = lower_level->add_rq(&VAPQ.front());
      VAPQ.pop_front();
    }
        if (result > 0)
        pf_issued++;

  } else {

    // TEMPORARY SOLUTION: mark prefetches as translated after a fixed latency
    if (VAPQ.has_ready()) {
      VAPQ.front().address = vmem.va_to_pa(cpu, VAPQ.front().v_address).first;

      // move the translated prefetch over to the regular PQ
      int result = add_pq(&VAPQ.front());

      // remove the prefetch from the VAPQ
      if (result != -2)
        VAPQ.pop_front();

      if (result > 0)
        pf_issued++;
    }
  }
}

int CACHE::add_pq(PACKET* packet)
{
  assert(packet->address != 0);
  PQ_ACCESS++;

  DP(if (warmup_complete[packet->cpu]) {
    std::cout << "[" << NAME << "_WQ] " << __func__ << " instr_id: " << packet->instr_id << " address: " << std::hex << (packet->address >> OFFSET_BITS);
    std::cout << " full_addr: " << packet->address << " v_address: " << packet->v_address << std::dec << " type: " << +packet->type
              << " occupancy: " << RQ.occupancy();
  })

  // check for the latest wirtebacks in the write queue
  champsim::delay_queue<PACKET>::iterator found_wq = std::find_if(WQ.begin(), WQ.end(), eq_addr<PACKET>(packet->address, match_offset_bits ? 0 : OFFSET_BITS));

  if (found_wq != WQ.end()) {

    DP(if (warmup_complete[packet->cpu]) std::cout << " MERGED_WQ" << std::endl;)

    packet->data = found_wq->data;
    for (auto ret : packet->to_return)
      ret->return_data(packet);

    WQ_FORWARD++;
    return -1;
  }

  // check for duplicates in the PQ
  auto found = std::find_if(PQ.begin(), PQ.end(), eq_addr<PACKET>(packet->address, OFFSET_BITS));

  if (found != PQ.end()) {
    DP(if (warmup_complete[packet->cpu]) std::cout << " MERGED_PQ" << std::endl;)

    found->fill_level = std::min(found->fill_level, packet->fill_level);
    packet_dep_merge(found->to_return, packet->to_return);

    PQ_MERGED++;
    return 0;
  }

  // check occupancy
  if (PQ.full()) {

    DP(if (warmup_complete[packet->cpu]) std::cout << " FULL" << std::endl;)

    PQ_FULL++;
    return -2; // cannot handle this request
  }

  // if there is no duplicate, add it to PQ
  if (warmup_complete[cpu]) {
    PQ.push_back(*packet);
  } else {
    PQ.push_back_ready(*packet);
  }

  // DP(if (warmup_complete[packet->cpu]) std::cout << " ADDED" << std::endl;)
  if(packet->is_stlb_prefetch == 1){
    // cout << "Added to PQ " << packet->address << endl;

  }
 

  PQ_TO_CACHE++;
  return (PQ.occupancy());
}

void CACHE::return_data(PACKET* packet)
{
  // check MSHR information
  auto mshr_entry = std::find_if(MSHR.begin(), MSHR.end(), eq_addr<PACKET>(packet->address, OFFSET_BITS));
  auto first_unreturned = std::find_if(MSHR.begin(), MSHR.end(), [](auto x) { return x.event_cycle == std::numeric_limits<uint64_t>::max(); });

  // sanity check
  if (mshr_entry == MSHR.end()) {
    std::cerr << "[" << NAME << "_MSHR] " << __func__ << " instr_id: " << packet->instr_id << " cannot find a matching entry!";
    std::cerr << " address: " << std::hex << packet->address;
    std::cerr << " v_address: " << packet->v_address;
    std::cerr << " address: " << (packet->address >> OFFSET_BITS) << std::dec;
    std::cerr << " event: " << packet->event_cycle << " current: " << current_cycle << std::endl;
    assert(0);
  }

  // MSHR holds the most updated information about this request
  mshr_entry->data = packet->data;
  mshr_entry->pf_metadata = packet->pf_metadata;
  mshr_entry->event_cycle = current_cycle + (warmup_complete[cpu] ? FILL_LATENCY : 0);

  DP(if (warmup_complete[packet->cpu]) {
    std::cout << "[" << NAME << "_MSHR] " << __func__ << " instr_id: " << mshr_entry->instr_id;
    std::cout << " address: " << std::hex << (mshr_entry->address >> OFFSET_BITS) << " full_addr: " << mshr_entry->address;
    std::cout << " data: " << mshr_entry->data << std::dec;
    std::cout << " index: " << std::distance(MSHR.begin(), mshr_entry) << " occupancy: " << get_occupancy(0, 0);
    std::cout << " event: " << mshr_entry->event_cycle << " current: " << current_cycle << std::endl;
  });

  // Order this entry after previously-returned entries, but before non-returned
  // entries
  std::iter_swap(mshr_entry, first_unreturned);
}

uint32_t CACHE::get_occupancy(uint8_t queue_type, uint64_t address)
{
  if (queue_type == 0)
    return std::count_if(MSHR.begin(), MSHR.end(), is_valid<PACKET>());
  else if (queue_type == 1)
    return RQ.occupancy();
  else if (queue_type == 2)
    return WQ.occupancy();
  else if (queue_type == 3)
    return PQ.occupancy();

  return 0;
}

uint32_t CACHE::get_size(uint8_t queue_type, uint64_t address)
{
  if (queue_type == 0)
    return MSHR_SIZE;
  else if (queue_type == 1)
    return RQ.size();
  else if (queue_type == 2)
    return WQ.size();
  else if (queue_type == 3)
    return PQ.size();

  return 0;
}

bool CACHE::should_activate_prefetcher(int type) { return (1 << static_cast<int>(type)) & pref_activate_mask; }

void CACHE::print_deadlock()
{
  if (!std::empty(MSHR)) {
    std::cout << NAME << " MSHR Entry" << std::endl;
    std::size_t j = 0;
    for (PACKET entry : MSHR) {
      std::cout << "[" << NAME << " MSHR] entry: " << j++ << " instr_id: " << entry.instr_id;
      std::cout << " address: " << std::hex << (entry.address >> LOG2_BLOCK_SIZE) << " full_addr: " << entry.address << std::dec << " type: " << +entry.type;
      std::cout << " fill_level: " << +entry.fill_level << " event_cycle: " << entry.event_cycle << std::endl;
    }
  } else {
    std::cout << NAME << " MSHR empty" << std::endl;
  }
}
