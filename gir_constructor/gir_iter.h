// Neighbour iteration over the current graph (base CSR minus tombstones plus
// per-vertex deltas). Shared by gir.cpp and gir_query.cpp.
#ifndef GIR_ITER_H
#define GIR_ITER_H
#include "gir.h"
#include "graph.h"

template <class F>
inline void GIR::forOut(int v, F&& f) const {
    const bool chk = n_dead_dir != 0;
    for (index_t i = g_directed->fw_beg_pos[v]; i < g_directed->fw_beg_pos[v + 1]; ++i) {
        int w = (int)g_directed->fw_csr[i];
        if (chk && bitGet(dead_fw, (long long)i)) continue;
        f(w);
    }
    for (int w : extra_out_adj[v]) f(w);
}

template <class F>
inline void GIR::forIn(int v, F&& f) const {
    const bool chk = n_dead_dir != 0;
    for (index_t i = g_directed->bw_beg_pos[v]; i < g_directed->bw_beg_pos[v + 1]; ++i) {
        int w = (int)g_directed->bw_csr[i];
        if (chk && bitGet(dead_bw, (long long)i)) continue;
        f(w);
    }
    for (int w : extra_in_adj[v]) f(w);
}

template <class F>
inline void GIR::forUnd(int v, F&& f) const {
    const bool chk = n_dead_und != 0;
    for (long long i = sym_beg[v]; i < sym_beg[v + 1]; ++i) {
        int w = sym_adj[i];
        if (chk && bitGet(dead_sym, i)) continue;
        f(w);
    }
    for (int w : extra_und_adj[v]) f(w);
}

template <class F>
inline bool GIR::anyUnd(int v, F&& f) const {
    const bool chk = n_dead_und != 0;
    for (long long i = sym_beg[v]; i < sym_beg[v + 1]; ++i) {
        if (chk && bitGet(dead_sym, i)) continue;
        if (f(sym_adj[i])) return true;
    }
    for (int w : extra_und_adj[v]) if (f(w)) return true;
    return false;
}

template <class F>
inline bool GIR::anyOut(int v, F&& f) const {
    const bool chk = n_dead_dir != 0;
    for (index_t i = g_directed->fw_beg_pos[v]; i < g_directed->fw_beg_pos[v + 1]; ++i) {
        if (chk && bitGet(dead_fw, (long long)i)) continue;
        if (f((int)g_directed->fw_csr[i])) return true;
    }
    for (int w : extra_out_adj[v]) if (f(w)) return true;
    return false;
}

template <class F>
inline bool GIR::anyIn(int v, F&& f) const {
    const bool chk = n_dead_dir != 0;
    for (index_t i = g_directed->bw_beg_pos[v]; i < g_directed->bw_beg_pos[v + 1]; ++i) {
        if (chk && bitGet(dead_bw, (long long)i)) continue;
        if (f((int)g_directed->bw_csr[i])) return true;
    }
    for (int w : extra_in_adj[v]) if (f(w)) return true;
    return false;
}

#endif  // GIR_ITER_H
