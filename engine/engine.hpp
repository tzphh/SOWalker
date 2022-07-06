#ifndef _GRAPH_ENGINE_H_
#define _GRAPH_ENGINE_H_

#include <functional>
#include "cache.hpp"
#include "schedule.hpp"
#include "util/timer.hpp"
#include "metrics/metrics.hpp"
#include "apps/secondorder.hpp"


class graph_engine {
public:
    graph_cache                         *cache;
    graph_walk                          *walk_manager;
    graph_driver                        *driver;
    graph_config                        *conf;
    graph_timer                         gtimer;
    std::vector<RandNum>                seeds;
    metrics &_m;

#ifdef PROF_STEPS
    size_t total_times;
    double sum_avg_steps;
#endif

    graph_engine(graph_cache& _cache, graph_walk& manager, graph_driver& _driver, graph_config& _conf, metrics &m) : _m(m){
        cache         = &_cache;
        walk_manager = &manager;
        driver        = &_driver;
        conf          = &_conf;
        seeds = std::vector<RandNum>(conf->nthreads, RandNum(9898676785859));
        for(tid_t tid = 0; tid < conf->nthreads; tid++) {
            seeds[tid] = time(NULL) + tid;
        }

#ifdef PROF_STEPS
        total_times = 0;
        sum_avg_steps = 0.0;
#endif

    }

    void prologue(second_order_app_t &userprogram, std::function<void(graph_walk *)> init_func = nullptr)
    {
        logstream(LOG_INFO) << "  =================  STARTED  ======================  " << std::endl;
        logstream(LOG_INFO) << "Random walks, random generate " << userprogram.get_numsources() << " walks on whole graph, exec_threads = " << conf->nthreads << std::endl;
        logstream(LOG_INFO) << "vertices : " << conf->nvertices << ", edges : " << conf->nedges << std::endl;
        srand(time(0));

        omp_set_num_threads(conf->nthreads);
        _m.start_time("run_app");
        userprogram.prologue(walk_manager, init_func);
    }

    void run(second_order_app_t &userprogram, scheduler *block_scheduler)
    {
        logstream(LOG_DEBUG) << "graph blocks : " << walk_manager->global_blocks->nblocks << ", memory blocks : " << cache->ncblock << std::endl;
        logstream(LOG_INFO) << "Random walks start executing, please wait for a minute." << std::endl;
        gtimer.start_time();
        int run_count = 0;
        wid_t interval_max_walks = conf->max_nthreads * MAX_TWALKS * 5;
        bid_t nblocks = walk_manager->nblocks;
        while(!walk_manager->test_finished_walks()) {
            wid_t total_walks = walk_manager->nwalks();
            logstream(LOG_DEBUG) << "run time : " << gtimer.runtime() << std::endl;
            logstream(LOG_DEBUG) << "run_count = " << run_count << ", total walks = " << total_walks << std::endl;
            cache->walk_blocks.clear();
            block_scheduler->schedule(*cache, *driver, *walk_manager);
            size_t pos = 0;
            logstream(LOG_DEBUG) << "cache walk block size : " << cache->walk_blocks.size() << std::endl;
            for(bid_t blk = 0; blk < walk_manager->global_blocks->nblocks; blk++) {
                std::cout << (*(walk_manager->global_blocks))[blk].cache_index << " ";
            }
            std::cout << std::endl;
            while(pos < cache->walk_blocks.size()) {
                wid_t nwalks = 0;
                walk_manager->walks.clear();
                while(pos < cache->walk_blocks.size() && nwalks + walk_manager->nmwalks(cache->walk_blocks[pos]) <= interval_max_walks) {
                    nwalks += walk_manager->nmwalks(cache->walk_blocks[pos]);
                    _m.start_time("walk I/O");
                    walk_manager->load_memory_walks(cache->walk_blocks[pos]);
                    _m.stop_time("walk I/O");
                    pos++;
                }
                update_walk(userprogram, nwalks);
                logstream(LOG_DEBUG) << "load memory walks, pos = " << pos << ", walks = " << nwalks << std::endl;
            }
            pos = 0;
            while (pos < cache->walk_blocks.size()) {
                wid_t num_disk_walks = walk_manager->ndwalks(cache->walk_blocks[pos]), disk_load_walks = 0;
                while(num_disk_walks > 0){
                    wid_t interval_walks = std::min(num_disk_walks, interval_max_walks);
                    num_disk_walks -= interval_walks;
                    _m.start_time("walk I/O");
                    wid_t nwalks = walk_manager->load_disk_walks(cache->walk_blocks[pos], interval_walks, disk_load_walks);
                    _m.stop_time("walk I/O");
                    disk_load_walks += interval_walks;
                    update_walk(userprogram, nwalks);
                    logstream(LOG_DEBUG) << "load disk walks from " << cache->walk_blocks[pos] / nblocks << " to " << cache->walk_blocks[pos] % nblocks << ", walks = " << nwalks << std::endl;
                }
                _m.start_time("walk I/O");
                walk_manager->dump_walks(cache->walk_blocks[pos]);
                _m.stop_time("walk I/O");
                pos++;
            }
            run_count++;
        }
        logstream(LOG_DEBUG) << gtimer.runtime() << "s, total run count : " << run_count << std::endl;
    }

    void epilogue(second_order_app_t &userprogram)
    {
        userprogram.epilogue();
        _m.stop_time("run_app");
#ifdef PROF_STEPS
        std::cout << "each walk step : " << sum_avg_steps / total_times << std::endl;
#endif
        logstream(LOG_INFO) << "  ================= FINISHED ======================  " << std::endl;
    }

    void update_walk(second_order_app_t &userprogram, wid_t nwalks)
    {
        if(nwalks < 100) omp_set_num_threads(1);
        else omp_set_num_threads(conf->nthreads);
        _m.start_time("exec_block_walk");
        {           
        #pragma omp parallel for schedule(dynamic)
            for(wid_t idx = 0; idx < nwalks; idx++) {
                userprogram.update_walk(walk_manager->walks[idx], cache, walk_manager, &seeds[omp_get_thread_num()]);
            }
#ifdef PROF_STEPS
            total_times++;
            sum_avg_steps += (double)run_steps / nwalks;        
#endif
        }
        _m.stop_time("exec_block_walk");
    }
};

#endif
