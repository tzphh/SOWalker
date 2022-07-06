#include <omp.h>
#include "api/constants.hpp"
#include "engine/config.hpp"
#include "engine/cache.hpp"
#include "engine/schedule.hpp"
#include "engine/walk.hpp"
#include "engine/engine.hpp"
#include "logger/logger.hpp"
#include "util/io.hpp"
#include "util/util.hpp"
#include "metrics/metrics.hpp"
#include "metrics/reporter.hpp"
#include "apps/secondorder.hpp"
#include "preprocess/graph_converter.hpp"

int main(int argc, const char *argv[])
{
    assert(argc >= 2);
    set_argc(argc, argv);
    logstream(LOG_INFO) << "app : " << argv[0] << ", dataset : " << argv[1] << std::endl;
    std::string input = argv[1];
    bool weighted = get_option_bool("weighted");
    bool sorted = get_option_bool("sorted");
    bool skip = get_option_bool("skip"); // use to skip the interactive convert query
    size_t blocksize = get_option_long("blocksize", BLOCK_SIZE);
    size_t nthreads = get_option_int("nthreads", omp_get_max_threads());
    size_t dynamic = get_option_bool("dynamic"); // the blocksize is dynamic, according to the number of walks
    size_t cache_size = get_option_int("cache", MEMORY_CACHE / (1024LL * 1024 * 1024));
    size_t max_iter = get_option_int("iter", 30);
    wid_t walks = (wid_t)get_option_int("walksource", 10000000);
    wid_t walkpersource = (wid_t)get_option_int("walkpersource", 1);
    hid_t steps = (hid_t)get_option_int("length", 20);
    real_t alpha = (real_t)get_option_float("alpha", 0.2);

    auto static_query_blocksize = [blocksize](vid_t nvertices) { return blocksize; };
    auto dynamic_query_blocksize = [&blocksize, walks](vid_t nvertices) {
        wid_t nwalks = 1LL * walks * nvertices;
        int dom = 0;
        while(nwalks){
            dom++;
            nwalks /= 10;
        }

        blocksize= 1LL * pow(2, 11 + dom) * 1024;
        return blocksize;
    };

    std::function<size_t(vid_t nvertices)> query_blocksize;
    if(dynamic) query_blocksize = dynamic_query_blocksize;
    else query_blocksize = static_query_blocksize;
    graph_converter converter(remove_extension(argv[1]), weighted, sorted);
    convert(input, converter, query_blocksize, skip);
    std::string base_name = remove_extension(argv[1]);

    /* graph meta info */
    vid_t nvertices;
    eid_t nedges;
    load_graph_meta(base_name, &nvertices, &nedges, weighted);

    graph_config conf = {
        base_name,
        cache_size * 1024LL * 1024 * 1024,
        blocksize,
        (tid_t)nthreads,
        (tid_t)omp_get_max_threads(),
        nvertices,
        nedges,
        weighted
    };
    
    graph_block blocks(&conf);
    bid_t nmblocks = get_option_int("nmblocks", blocks.nblocks);

    metrics m("autoregressive_walksource_" + std::to_string(walks) + "_steps_" + std::to_string(steps) + "_dataset_" + argv[1] + "_iter_" + std::to_string(max_iter));
    graph_driver driver(&conf, m);

    graph_walk walk_mangager(conf, driver, blocks);
    graph_cache cache(min_value(nmblocks, blocks.nblocks), &conf);
    m.set("nblocks", std::to_string(blocks.nblocks));
    m.set("ncblocks", std::to_string(cache.ncblock));
    
    autoregreesive_app_t userprogram(walks, steps, false,alpha);
    //autoregreesive_app_t userprogram(walks, steps, alpha);
    graph_engine engine(cache, walk_mangager, driver, conf, m);

    //lp_solver_scheduler_t walk_scheduler(m);
    //simulated_annealing_scheduler_t walk_scheduler(max_iter, m);
    navie_graphwalker_scheduler_t walk_scheduler(m);
    // random_scheduler_t walk_scheduler(m);
    // greedy_scheduler_t walk_scheduler(m);
    // greedy_graphwalker_scheduler_t walk_scheduler(m);

    auto init_func = [walks, walkpersource](graph_walk *walk_manager)
    {
        #pragma omp parallel for schedule(static)
        for(wid_t off = 0; off < walks; off++)
        {
            vid_t vertex = rand() % walk_manager->nvertices;
            bid_t index = walk_manager->global_blocks->get_block(vertex);
            for(wid_t cnt = 0; cnt < walkpersource; cnt++) {
                wid_t walker_id = off * walkpersource + cnt;
                walker_t walker = walker_makeup(walker_id, vertex, vertex, vertex, 0, index, index);
                walk_manager->move_walk(walker);
              }
        }
    };
    
    engine.prologue(userprogram, init_func);
    engine.run(userprogram, &walk_scheduler);
    engine.epilogue(userprogram);

    metrics_report(m);
    return 0;
}
