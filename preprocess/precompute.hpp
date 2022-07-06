#ifndef _GRAPH_PRECOMPUTE_H_
#define _GRAPH_PRECOMPUTE_H_

#include <string>
#include <vector>
#include <numeric>
#include <queue>
#include <omp.h>
#include <algorithm>
#include <unordered_set>
#include "api/types.hpp"
#include "util/util.hpp"
#include "util/io.hpp"
#include "util/hash.hpp"

/* the block structure used for precompute */
struct pre_block_t {
    vid_t start_vert, nverts;
    eid_t start_edge, nedges;
    eid_t *beg_pos;
    vid_t *csr;
    real_t *weights;

    pre_block_t() {
        beg_pos = NULL;
        csr = NULL;
        weights = NULL;
    }
    ~pre_block_t() {
        if(beg_pos) free(beg_pos);
        if(csr) free(csr);
        if(weights) free(weights);
    }
};


static void sort_block_vertex_neighbors(const pre_block_t& block, vid_t *new_csr, real_t *new_weights)
{
    vid_t *csr_start = block.csr;
    auto cmp = [&csr_start](eid_t u, eid_t v) { return csr_start[u] < csr_start[v]; };
    omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel for schedule(static)
    for(vid_t vertex = 0; vertex < block.nverts; ++vertex) {
        eid_t adj_head = block.beg_pos[vertex] - block.start_edge, adj_tail = block.beg_pos[vertex + 1] - block.start_edge;
        std::vector<eid_t> index(adj_tail - adj_head);
        std::iota(index.begin(), index.end(), adj_head);
        std::sort(index.begin(), index.end(), cmp);
        for(eid_t off = adj_head; off < adj_tail; ++off) {
            new_csr[off] = block.csr[index[off - adj_head]];
            if(new_weights != nullptr) {
                new_weights[off] = block.weights[index[off - adj_head]];
            }
        }
    }
}

void sort_vertex_neighbors(const std::string &base_name, size_t blocksize, bool weighted)
{
    std::string vert_block_name = get_vert_blocks_name(base_name, blocksize);
    std::string edge_block_name = get_edge_blocks_name(base_name, blocksize);
    std::string weight_name = get_weights_name(base_name);
    std::string beg_pos_name = get_beg_pos_name(base_name);
    std::string csr_name = get_csr_name(base_name);

    std::vector<vid_t> vblocks = load_graph_blocks<vid_t>(vert_block_name);
    std::vector<eid_t> eblocks = load_graph_blocks<eid_t>(edge_block_name);

    int vertdesc = open(beg_pos_name.c_str(), O_RDONLY);
    int edgedesc = open(csr_name.c_str(), O_RDWR);
    int whtdesc = 0;
    if(weighted) {
        whtdesc = open(weight_name.c_str(), O_RDWR);
    }
    assert(vertdesc > 0 && edgedesc > 0);
    bid_t nblocks = vblocks.size() - 1;
    logstream(LOG_INFO) << "load vblocks and eblocks successfully, block count : " << nblocks << std::endl;
    pre_block_t block;
    logstream(LOG_INFO) << "start to sort the vertex neighbors, nblocks = " << nblocks << std::endl;

    vid_t *new_csr = nullptr;
    real_t *new_weights = nullptr;
    for (bid_t blk = 0; blk < nblocks; blk++)
    {
        block.nverts = vblocks[blk + 1] - vblocks[blk];
        block.nedges = eblocks[blk + 1] - eblocks[blk];
        block.start_vert = vblocks[blk];
        block.start_edge = eblocks[blk];

        block.beg_pos = (eid_t *)realloc(block.beg_pos, (block.nverts + 1) * sizeof(eid_t));
        block.csr = (vid_t *)realloc(block.csr, block.nedges * sizeof(vid_t));
        if(weighted) {
            block.weights = (real_t *)realloc(block.weights, block.nedges * sizeof(real_t));
        }

        load_block_range(vertdesc, block.beg_pos, (block.nverts + 1), block.start_vert * sizeof(eid_t));
        load_block_range(edgedesc, block.csr, block.nedges, block.start_edge * sizeof(vid_t));
        if(weighted) {
            load_block_range(whtdesc, block.weights, block.nedges, block.start_edge * sizeof(real_t));
        }

        new_csr = (vid_t *)realloc(new_csr, block.nedges * sizeof(vid_t));
        if(weighted) {
            new_weights = (real_t *)realloc(new_weights, block.nedges * sizeof(real_t));
        }

        sort_block_vertex_neighbors(block, new_csr, new_weights);
        dump_block_range(edgedesc, new_csr, block.nedges, block.start_edge * sizeof(vid_t));
        if(weighted) {
            dump_block_range(whtdesc, new_weights, block.nedges, block.start_edge * sizeof(real_t));
        }
        logstream(LOG_INFO) << "finish sort vertex neighbors for block = " << blk << std::endl;
    }

    if(new_csr) free(new_csr);
    if(new_weights) free(new_weights);
}

real_t calc_block_expect_walk_len(pre_block_t *block, size_t len_limit)
{
    std::unordered_set<vid_t> tmp_block_vertices;
    std::vector<vid_t> block_vertics(block->nverts);
    std::iota(block_vertics.begin(), block_vertics.end(), block->start_vert);
    std::vector<real_t> access_wht(block->nverts, 1.0 / block->nverts), bak_access_wht(block->nverts, 0.0);
    size_t len = 1;

    real_t exp_walk_len = 1.0, base = 1.0;

    while (len <= len_limit)
    {
        real_t r = 0.0;
        for (auto v : block_vertics)
        {
            vid_t v_off = v - block->start_vert;
            eid_t adj_start = block->beg_pos[v_off] - block->start_edge, adj_end = block->beg_pos[v_off + 1] - block->start_edge;
            size_t inner_count = 0;
            for (eid_t off = adj_start; off < adj_end; off++)
            {
                if (block->csr[off] >= block->start_vert && block->csr[off] < block->start_vert + block->nverts)
                {
                    tmp_block_vertices.insert(block->csr[off]);
                    inner_count++;
                    bak_access_wht[block->csr[off] - block->start_vert] += access_wht[v_off] / (adj_end - adj_start);
                }
            }

            if (adj_end > adj_start && inner_count > 0)
            {
                r += access_wht[v_off] * inner_count / (adj_end - adj_start);
            }
        }
        base *= r;
        exp_walk_len += len * base;
        logstream(LOG_DEBUG) << "len :  " << len << ", next step in block prob :  " << r << ", joint prob : " << base << std::endl;
        logstream(LOG_DEBUG) << "cal block expect walk len : " << exp_walk_len << ", tmp_block_vertices size : " << tmp_block_vertices.size() << std::endl;
        len++;

        if (tmp_block_vertices.empty())
            break;
        block_vertics = std::vector<vid_t>(tmp_block_vertices.begin(), tmp_block_vertices.end());
        tmp_block_vertices.clear();

        /* calculate the access weight for each vertex */
        real_t wht_sum = std::accumulate(bak_access_wht.begin(), bak_access_wht.end(), 0.0);
        for (auto &wht : bak_access_wht)
            wht /= wht_sum;
        access_wht = bak_access_wht;
        std::fill(bak_access_wht.begin(), bak_access_wht.end(), 0.0);
    }
    return exp_walk_len;
}

/**
 * This method does following thing:
 * compute the expected walk length for each block
 */
void calc_expected_walk_length(const std::string &base_name, size_t blocksize, size_t len_limit)
{
    std::string vert_block_name = get_vert_blocks_name(base_name, blocksize);
    std::string edge_block_name = get_edge_blocks_name(base_name, blocksize);
    std::string csr_name = get_csr_name(base_name);
    std::string beg_pos_name = get_beg_pos_name(base_name);

    std::vector<vid_t> vblocks = load_graph_blocks<vid_t>(vert_block_name);
    std::vector<eid_t> eblocks = load_graph_blocks<eid_t>(edge_block_name);

    int vertdesc = open(beg_pos_name.c_str(), O_RDONLY);
    int edgedesc = open(csr_name.c_str(), O_RDONLY);

    bid_t nblocks = vblocks.size() - 1;
    logstream(LOG_INFO) << "load vblocks and eblocks successfully, block count : " << nblocks << std::endl;
    pre_block_t block;
    logstream(LOG_INFO) << "start to compute block expected walk length, nblocks = " << nblocks << std::endl;
    std::vector<real_t> block_walk_len(nblocks, 0.0);

//#pragma omp parallel for schedule(dynamic)
    for (bid_t blk = 0; blk < nblocks; blk++)
    {
        block.nverts = vblocks[blk + 1] - vblocks[blk];
        block.nedges = eblocks[blk + 1] - eblocks[blk];
        block.start_edge = eblocks[blk];
        block.start_vert = vblocks[blk];
        block.beg_pos = (eid_t*)realloc(block.beg_pos, (block.nverts + 1) * sizeof(eid_t));
        load_block_range(vertdesc, block.beg_pos, block.nverts + 1, block.start_vert * sizeof(eid_t));
        block.csr = (vid_t *)realloc(block.csr, block.nedges * sizeof(vid_t));
        load_block_range(edgedesc, block.csr, block.nedges, block.start_edge * sizeof(vid_t));
        logstream(LOG_INFO) << "start computing expected walk length for block = " << blk << std::endl;     
        //#pragma omp critical 
            block_walk_len[blk] = calc_block_expect_walk_len(&block, len_limit);
        logstream(LOG_INFO) << "finish computing expected walk length for block = " << blk << std::endl;
    }
    close(edgedesc);
    std::string walk_length_name = get_expected_walk_length_name(base_name, blocksize);
    auto stream = std::fstream(walk_length_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    stream.write(reinterpret_cast<char*>(block_walk_len.data()), sizeof(real_t) * block_walk_len.size());
    stream.close();
}

#endif
