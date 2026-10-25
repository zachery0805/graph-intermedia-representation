// Converts a GIR protocol directory into GraphBolt/DZiG inputs:
//   loaded.adj / loaded_w.adj   (Weighted)AdjacencyGraph of the 50% graph, n = V
//   {insert,delete,mixed}.stream / *_w.stream   "a|d <TAB> u <TAB> v [<TAB> w]"
// Weights are GIR's edgeWeight (proto::weight). Usage: prepare <protocol dir> <out dir>
#include <algorithm>
#include <string>
#include "../common/protocol.h"

int main(int argc, char** argv) {
    if (argc < 3) return 2;
    const std::string dir = argv[1], out = argv[2];
    const uint32_t n = (uint32_t)proto::meta(dir, "V");
    auto el = proto::read_el(dir + "/loaded.el");
    std::sort(el.begin(), el.end());
    el.erase(std::unique(el.begin(), el.end()), el.end());
    std::vector<uint64_t> off(n + 1, 0);
    for (auto& e : el) off[e.first + 1]++;
    for (uint32_t i = 0; i < n; ++i) off[i + 1] += off[i];
    for (int w = 0; w < 2; ++w) {
        FILE* f = std::fopen((out + (w ? "/loaded_w.adj" : "/loaded.adj")).c_str(), "w");
        std::fprintf(f, "%s\n%u\n%zu\n", w ? "WeightedAdjacencyGraph" : "AdjacencyGraph", n, el.size());
        for (uint32_t i = 0; i < n; ++i) std::fprintf(f, "%lu\n", off[i]);
        for (auto& e : el) std::fprintf(f, "%u\n", e.second);
        if (w) for (auto& e : el) std::fprintf(f, "%u\n", proto::weight(e.first, e.second));
        std::fclose(f);
    }
    for (const char* s : {"insert", "delete", "mixed"}) {
        auto ops = proto::read_ops(dir + "/" + s + ".ops");
        for (int w = 0; w < 2; ++w) {
            FILE* f = std::fopen((out + "/" + s + (w ? "_w" : "") + ".stream").c_str(), "w");
            for (auto& o : ops) {
                if (w) std::fprintf(f, "%c\t%u\t%u\t%u\n", o.op, o.u, o.v, proto::weight(o.u, o.v));
                else std::fprintf(f, "%c\t%u\t%u\n", o.op, o.u, o.v);
            }
            std::fclose(f);
        }
    }
    std::printf("n=%u m=%zu\n", n, el.size());
    return 0;
}
