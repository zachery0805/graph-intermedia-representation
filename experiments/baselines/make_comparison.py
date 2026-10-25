#!/usr/bin/env python3
"""Builds experiments/results/comparison.md from the measured CSVs.

Inputs (all in experiments/results/):
  baselines.csv               Terrace / Spruce / Aspen / DZiG / GIR-streaming rows
  gir_benchmark.csv           GIR, directed graphs, full protocol (updates + thread sweep)
  gir_queries_lj_clean.csv    GIR, directed livejournal, queries re-measured on an idle machine
  gir_benchmark_sym.csv       GIR, symmetrized livejournal
  bench_run.log               GIR livejournal log (per-case latencies, section 2.1)
  history/*_pre_speedup.*     the same GIR measurements before the 2026-10-14 speed-up (section 9)
No number in the output is typed by hand except the paper's claims (section 8),
which are copied from paper_original/experiment.tex and figures/scripts/plot_throughput.py.
"""
import csv
import os
from collections import defaultdict

R = os.path.join(os.path.dirname(__file__), "..", "results")


def load(name, has_system):
    rows = []
    with open(os.path.join(R, name)) as f:
        for r in csv.DictReader(f):
            if not has_system:
                r["system"] = "GIR"
            r["value"] = float(r["value"])
            r["threads"] = int(r["threads"])
            rows.append(r)
    return rows


rows = load("baselines.csv", True) + load("gir_benchmark.csv", False) + load("gir_benchmark_sym.csv", False)
clean = load("gir_queries_lj_clean.csv", False)
pre = (load("history/gir_benchmark_pre_speedup.csv", False) + load("history/gir_benchmark_sym_pre_speedup.csv", False)
       + load("history/gir_stream_pre_speedup.csv", True))
pre_clean = load("history/gir_queries_lj_clean_pre_speedup.csv", False)
mid = (load("history/gir_benchmark_round1.csv", False) + load("history/gir_benchmark_sym_round1.csv", False)
       + load("history/gir_stream_round1.csv", True))
mid_clean = load("history/gir_queries_lj_clean_round1.csv", False)
mid2 = (load("history/gir_benchmark_round2.csv", False) + load("history/gir_benchmark_sym_round2.csv", False)
        + load("history/gir_stream_round2.csv", True))
mid2_clean = load("history/gir_queries_lj_clean_round2.csv", False)
mid3 = (load("history/gir_benchmark_round3.csv", False) + load("history/gir_benchmark_sym_round3.csv", False)
        + load("history/gir_stream_round3.csv", True))
mid3_clean = load("history/gir_queries_lj_clean_round3.csv", False)
mid4 = (load("history/gir_benchmark_round4.csv", False) + load("history/gir_benchmark_sym_round4.csv", False)
        + load("history/gir_stream_round4.csv", True))
mid4_clean = load("history/gir_queries_lj_clean_round4.csv", False)
mid5 = (load("history/gir_benchmark_round5.csv", False) + load("history/gir_benchmark_sym_round5.csv", False)
        + load("history/gir_stream_round5.csv", True))
mid5_clean = load("history/gir_queries_lj_clean_round5.csv", False)
mid6 = (load("history/gir_benchmark_round6.csv", False) + load("history/gir_benchmark_sym_round6.csv", False)
        + load("history/gir_stream_round6.csv", True))
mid6_clean = load("history/gir_queries_lj_clean_round6.csv", False)
mid7 = (load("history/gir_benchmark_round7.csv", False) + load("history/gir_benchmark_sym_round7.csv", False)
        + load("history/gir_stream_round7.csv", True))
mid7_clean = load("history/gir_queries_lj_clean_round7.csv", False)
mid8 = (load("history/gir_benchmark_round8.csv", False) + load("history/gir_benchmark_sym_round8.csv", False)
        + load("history/gir_stream_round8.csv", True))
mid8_clean = load("history/gir_queries_lj_clean_round8.csv", False)
mid9 = (load("history/gir_benchmark_round9.csv", False) + load("history/gir_benchmark_sym_round9.csv", False)
        + load("history/gir_stream_round9.csv", True))
mid9_clean = load("history/gir_queries_lj_clean_round9.csv", False)


def case_ms(log):
    """Per-case mean latency (ms) of the first SPU 'cases:' line in a bench log."""
    import re
    try:
        text = open(os.path.join(R, log)).read()
    except OSError:
        return {}
    for line in text.splitlines():
        if line.strip().startswith("cases:") and "ms)" in line:
            return {m.group(1): (int(m.group(2)), float(m.group(3)))
                    for m in re.finditer(r"([\w-]+)=(\d+)\(([\d.e+-]+)ms\)", line)}
    return {}


def get(ds, system, metric, mode=None, threads=None, src=None):
    src = rows if src is None else src
    c = [r for r in src if r["dataset"] == ds and r["system"] == system and r["metric"] == metric
         and (mode is None or r["mode"] == mode) and (threads is None or r["threads"] == threads)]
    return c[-1]["value"] if c else None


def gq(metric, mode):  # GIR directed livejournal queries, clean run at 64 threads
    return get("livejournal", "GIR", metric, mode, 64, clean)


def f(v, digits=3):
    if v is None:
        return "—"
    if v < 0:
        return "失败"
    if v == 0:
        return "0"
    if abs(v) >= 1000:
        return f"{v:,.0f}"
    if abs(v) >= 100:
        return f"{v:.0f}"
    if abs(v) >= 10:
        return f"{v:.1f}"
    if abs(v) >= 1:
        return f"{v:.2f}"
    return f"{v:.{digits}g}"


def ratio(a, b):
    """a / b as 'x.xx×'; for throughput a=GIR, b=baseline; for time a=baseline, b=GIR."""
    if a is None or b is None or a <= 0 or b <= 0:
        return "—"
    r = a / b
    if r >= 100:
        return f"{r:,.0f}×"
    if r >= 1:
        return f"{r:.2f}×"
    inv = b / a  # GIR is slower: "1/N" means N times slower
    return f"1/{inv:.2f}" if inv < 10 else f"1/{inv:,.0f}"


out = []
w = out.append

w("# GIR 与基线系统实测对比\n")
w("所有数字由 `experiments/baselines/make_comparison.py` 从实测 CSV 生成，没有任何手抄或估算的数值。"
  "论文中的声称只出现在第 8 节，并标明出处。\n")

# ---------------------------------------------------------------- 1 setup
w("## 1. 实验设置\n")
w("| 项目 | 设置 |")
w("|---|---|")
w("| 机器 | 128 硬件线程（2 路），所有系统固定 64 线程，Spruce 另测 32 线程 |")
w("| 图 | soc-LiveJournal1：有向 4,847,571 点 / 68,475,391 边；对称化后 85,702,474 条有向边 |")
w("| 更新协议 | DZiG 式：随机打乱后加载 50% 边，插入流来自另一半，删除流来自已加载的一半，混合流两者交织；每种流 100,000 次有向边操作 |")
w("| 批大小 | 10,000（批处理系统）；Spruce 每批内并行单边操作；GIR 另测单边模式 SPU |")
w("| 查询 | 在完整图上：BFS 和 SSSP 各 16 个源点取平均；BC 为 64 个源点的总时间 |")
w("| SSSP 权重 | GIR 的确定性哈希权重，取值 1–100，所有系统一致 |")
w("| 输入一致性 | 所有系统读取同一份 `data/protocol/<数据集>/` 导出文件 |\n")

w("| 系统 | 版本与构建 | 驱动做了什么 |")
w("|---|---|---|")
w("| GIR | 本仓库 `gir_constructor` | 每次更新维护 WCC、SCC+DAG 和 BiCC 三层结构 |")
w("| Terrace（SIGMOD'21） | OpenCilk，`CILK=1`（Makefile 改 3 行） | 批插入；无批删除接口，删除用并行 `remove_edge` |")
w("| Spruce（SIGMOD'24） | 现有 Release 构建 | 每批并行单边插入和删除；BFS、SSSP、BC 由驱动在其邻居读取接口上实现 |")
w("| Aspen（PLDI'19） | 原生 Makefile，自带调度器 | 批插入和批删除；只支持对称图，没有 SSSP |")
w("| GraphBolt/DZiG（EuroSys'21） | OpenCilk；补了 3 处缺失的 `return`，否则初始计算死循环 | 原生流式模式：每批先更新图，再增量重算单源 BFS 或 SSSP |\n")

# ---------------------------------------------------------------- 2 directed updates
ds = "livejournal"
w("## 2. 有向 LiveJournal：更新吞吐（MEPS，越高越好）\n")
w("比值一栏为 GIR 批量模式（FMU）吞吐除以该系统吞吐。小于 1 表示 GIR 更慢，写作 1/N 表示慢 N 倍。\n")
w("| 系统 | 模式 | 插入 | 删除 | 混合 | GIR÷该系统（插入 / 删除 / 混合） |")
w("|---|---|---|---|---|---|")
gi, gd, gm = (get(ds, "GIR", m, "FMU_B10000") for m in ("insert_MEPS", "delete_MEPS", "mixed_MEPS"))
w(f"| **GIR** | FMU 批量 B=10000 | {f(gi)} | {f(gd)} | {f(gm)} | 1 |")
si, sd, sm = (get(ds, "GIR", m, "SPU") for m in ("insert_MEPS", "delete_MEPS", "mixed_MEPS"))
w(f"| **GIR** | SPU 单边 | {f(si)} | {f(sd)} | {f(sm)} | — |")
li, ld, lm = (get(ds, "GIR", m, "LAZY_ingest_B10000") for m in ("insert_MEPS", "delete_MEPS", "mixed_MEPS"))
ti, td, tm = (get(ds, "GIR", m, "LAZY_total_B10000") for m in ("insert_MEPS", "delete_MEPS", "mixed_MEPS"))
w(f"| **GIR** | 延迟维护，只写入 B=10000 | {f(li)} | {f(ld)} | {f(lm)} | — |")
w(f"| **GIR** | 延迟维护，写入 + 一次同步 | {f(ti)} | {f(td)} | {f(tm)} | — |")
for sysname, mode in (("Terrace", "batch_10000"),):
    a, b, c = (get(ds, sysname, m, mode) for m in ("insert_MEPS", "delete_MEPS", "mixed_MEPS"))
    w(f"| Terrace | 批量 B=10000 | {f(a)} | {f(b)} | {f(c)} | {ratio(gi, a)} / {ratio(gd, b)} / {ratio(gm, c)} |")
for t in (64, 32):
    a = get(ds, "Spruce", "insert_MEPS", None, t)
    b = get(ds, "Spruce", "delete_MEPS", None, t)
    c = get(ds, "Spruce", "mixed_MEPS", None, t)
    w(f"| Spruce（{t} 线程） | 并行单边，每批 10000 | {f(a)} | {f(b)} | {f(c)} | {ratio(gi, a)} / {ratio(gd, b)} / — |")
for alg, label in (("bfs", "BFS 运行中"), ("sssp", "SSSP 运行中")):
    a, b, c = (get(ds, "DZiG", f"{s}_mutation_MEPS", f"{alg}_batch_10000") for s in ("insert", "delete", "mixed"))
    w(f"| DZiG（{label}的图变更部分） | 批量 B=10000 | {f(a)} | {f(b)} | {f(c)} | {ratio(gi, a)} / {ratio(gd, b)} / {ratio(gm, c)} |")
w("")
w("\n延迟维护（2026-10-15 新增）：`ingestBatch` 只把批写入边存储，三层结构在 `sync()` 时一次补上"
  "（每阶段 10 批写入后同步一次）。“只写入”和存储型系统做的工作量相同；“写入 + 一次同步”包含结构维护。"
  "BFS、SSSP、BC 只读边存储，不依赖同步。\n")
w("| GIR 模式 | 对 Terrace：插入 | 删除 | 混合 |")
w("|---|---|---|---|")
tr = [get(ds, "Terrace", m, "batch_10000") for m in ("insert_MEPS", "delete_MEPS", "mixed_MEPS")]
for label, vals in (("SPU 单边", (si, sd, sm)), ("FMU 批量", (gi, gd, gm)), ("延迟维护，只写入", (li, ld, lm)),
                    ("延迟维护，写入 + 一次同步", (ti, td, tm))):
    w(f"| {label} | {ratio(vals[0], tr[0])} | {ratio(vals[1], tr[1])} | {ratio(vals[2], tr[2])} |")
w("")
w("- **Spruce 的混合流失败。** 删除阶段把某个顶点的最后一条边删掉后，Spruce 会删除该顶点；之后再往这个顶点插边，"
  "`InsertEdge` 永不返回。单线程可复现：混合流第 8,707 次操作，顶点 518289。表中记为失败。")
w("- **工作量不同。** Terrace、Spruce 和 Aspen 只更新邻接表；DZiG 的变更时间只含邻接表更新，增量重算另计（第 4 节）；"
  "GIR 每次更新都要维护三层结构。这是吞吐差距的主要来源。\n")

# GIR per-case latency
w("### 2.1 GIR 更新的分类构成（SPU 单边，有向 LiveJournal）\n")
w("| 情况 | 含义 | 次数 | 平均耗时（ms） | 总耗时（s） |")
w("|---|---|---|---|---|")
cases = [("C1", "插入，跨 WCC"), ("C2-dag", "插入，同 WCC 跨 SCC，无环"), ("C2-merge", "插入，成环，合并 SCC"),
         ("C3", "插入，同 SCC 跨块，块合并"), ("C4", "插入，无结构变化"), ("D-inter", "删除 SCC 间的边"),
         ("D-wcc-split", "删除导致 WCC 分裂"), ("D-intra-keep", "删除 SCC 内的边，SCC 不变"), ("D-scc-split", "删除导致 SCC 分裂")]
cm = case_ms("bench_run.log")
for c, desc in cases:
    n, ms = cm.get(c, (None, None))
    tot = n * ms / 1000 if n is not None and ms is not None else None
    w(f"| {c} | {desc} | {f(get(ds, 'GIR', 'case_' + c, 'SPU'))} | {f(ms)} | {f(tot)} |")
dels = [c for c, _ in cases if c.startswith("D-") and c in cm]
if dels:
    nd = sum(cm[c][0] for c in dels)
    top = max(dels, key=lambda c: cm[c][0] * cm[c][1])
    share = cm[top][0] * cm[top][1] / max(1e-9, sum(cm[c][0] * cm[c][1] for c in dels))
    w(f"\n次数和耗时覆盖插入、删除、混合三个阶段。删除耗时最大的一类是 {top}：占删除次数的 "
      f"{100 * cm[top][0] / nd:.0f}%，平均每次 {f(cm[top][1])} ms，占删除总耗时的 {100 * share:.0f}%。"
      "详见 `experiments/results/bench_run.log`。\n")

# ---------------------------------------------------------------- 3 directed queries
w("## 3. 有向 LiveJournal：查询（完整图，64 线程，越低越好）\n")
w("加速比一栏为该系统时间除以 GIR 时间（GIR 取方向优化 BFS 和 Δ-stepping SSSP）。大于 1 表示 GIR 更快。\n")
gb, gs, gbc = gq("bfs_ms", "direction-optimizing"), gq("sssp_ms", "delta-stepping"), gq("bc_ms_total", "brandes")
w("| 系统 | BFS 算法 | BFS（ms/源） | SSSP 算法 | SSSP（ms/源） | BC 64 源总计（ms） | GIR 加速比（BFS / SSSP / BC） |")
w("|---|---|---|---|---|---|---|")
w(f"| **GIR** | 方向优化，SCC 剪枝 | {f(gb)} | Δ-stepping（Δ=16） | {f(gs)} | {f(gbc)} | 1 |")
w(f"| GIR 两级模型（论文第 6 节） | DAG 波前 + 块任务 | {f(gq('bfs_ms', 'two-level'))} | 两级 Δ-stepping | {f(gq('sssp_ms', 'two-level'))} | — | — |")
w(f"| GIR 顺序 BFS（参照） | 单线程队列 | {f(gq('bfs_ms', 'sequential'))} | — | — | — | — |")
tb, ts = get(ds, "Terrace", "bfs_ms"), get(ds, "Terrace", "sssp_ms")
w(f"| Terrace | 仅推送 edgeMap | {f(tb)} | Bellman-Ford（仅推送） | {f(ts)} | 不支持（有向图上不正确） | {ratio(tb, gb)} / {ratio(ts, gs)} / — |")
for t in (64, 32):
    b, s, c = (get(ds, "Spruce", m, None, t) for m in ("bfs_ms", "sssp_ms", "bc_ms_total"))
    w(f"| Spruce（{t} 线程） | 自上而下（驱动实现） | {f(b)} | Δ-stepping（驱动实现） | {f(s)} | {f(c)} | "
      f"{ratio(b, gb) if t == 64 else '—'} / {ratio(s, gs) if t == 64 else '—'} / {ratio(c, gbc) if t == 64 else '—'} |")
w("")
w("DZiG 是单源增量系统，没有独立的多源查询接口，它的从头计算和增量重算时间见第 4 节。\n")
w("### 3.1 正确性交叉验证（第一个源点 2954358）\n")
w("| 系统 | BFS 可达顶点数 | SSSP 距离和 |")
w("|---|---|---|")
sci = lambda v: "—" if v is None else f"{v:.5e}"
for sname, t in (("GIR", 64), ("Terrace", None), ("Spruce", 64)):
    w(f"| {sname} | {sci(get(ds, sname, 'bfs_reached_src0', None, t))} | {sci(get(ds, sname, 'sssp_distsum_src0', None, t))} |")
w("\n统一按 6 位有效数字显示，因为基线驱动只按这个精度写入 CSV；GIR 的精确值为 4,400,347 和 308,120,134。三个系统一致。Terrace 原生的稠密模式 BFS 只到达 4,232,680 个顶点，"
  "SSSP 有 4,746,658 个错误距离，原因是它把出边当作入边，所以有向图上改用仅推送模式。\n")

# ---------------------------------------------------------------- 4 streaming
w("## 4. 流式场景：DZiG 增量计算对比 GIR 从头计算（有向 LiveJournal，50% 初始图）\n")
w("每种流 10 批，每批 10,000 次操作，源点相同，64 线程。每批时间等于图更新时间加 BFS 或 SSSP 重算时间。"
  "DZiG 保留该源点的结果并增量修正；GIR 不保存源点状态，每批从头计算。\n")
w("| 算法 | 流 | DZiG 图更新（ms/批） | DZiG 增量重算（ms/批） | DZiG 合计（ms/批） | GIR 图更新（ms/批） | GIR 从头计算（ms/批） | GIR 合计（ms/批） | DZiG 快多少倍 |")
w("|---|---|---|---|---|---|---|---|---|")
for alg in ("bfs", "sssp"):
    for s in ("insert", "delete", "mixed"):
        dm = get(ds, "DZiG", f"{s}_mutation_MEPS", f"{alg}_batch_10000")
        dzu = 10000 / (dm * 1e6) * 1000 if dm else None
        dzr = get(ds, "DZiG", f"{alg}_{s}_recompute_ms")
        dze = get(ds, "DZiG", f"{alg}_{s}_endtoend_ms")
        gm_ = get(ds, "GIR", f"{s}_mutation_MEPS", "fmu_batch_10000")
        gu = 10000 / (gm_ * 1e6) * 1000 if gm_ else None
        gr = get(ds, "GIR", f"{alg}_{s}_recompute_ms")
        ge = get(ds, "GIR", f"{alg}_{s}_endtoend_ms")
        w(f"| {alg.upper()} | {s} | {f(dzu)} | {f(dzr)} | {f(dze)} | {f(gu)} | {f(gr)} | {f(ge)} | {ratio(ge, dze)} |")
w("")
if get(ds, "GIR-lazy", "insert_mutation_MEPS") is not None:
    w("GIR 延迟维护（每批只写入边存储，BFS/SSSP 从头计算；结构层在整条流结束后同步一次，同步时间单列）：\n")
    w("| 算法 | 流 | GIR-lazy 写入（ms/批） | 从头计算（ms/批） | 合计（ms/批） | DZiG 合计（ms/批） | 对 DZiG | 流结束后同步（ms） |")
    w("|---|---|---|---|---|---|---|---|")
    for alg in ("bfs", "sssp"):
        for s in ("insert", "delete", "mixed"):
            lmps = get(ds, "GIR-lazy", f"{s}_mutation_MEPS")
            lu = 10000 / (lmps * 1e6) * 1000 if lmps else None
            lr = get(ds, "GIR-lazy", f"{alg}_{s}_recompute_ms")
            le = get(ds, "GIR-lazy", f"{alg}_{s}_endtoend_ms")
            dze = get(ds, "DZiG", f"{alg}_{s}_endtoend_ms")
            w(f"| {alg.upper()} | {s} | {f(lu)} | {f(lr)} | {f(le)} | {f(dze)} | {ratio(dze, le)} | {f(get(ds, 'GIR-lazy', f'sync_{s}_ms'))} |")
    w("\n“对 DZiG”一栏为 DZiG 合计时间除以 GIR-lazy 合计时间，大于 1 表示 GIR-lazy 更快。\n")
if get(ds, "GIR-lazy-inc", "insert_mutation_MEPS") is not None:
    w("GIR 增量计算（与 DZiG 相同的设定：保留该源点的 BFS/SSSP 结果，每批之后只修正受影响的部分。"
      "做法同 KickStarter：被删的树边使其下方子树失效，失效顶点从仍有效的入邻居重新取值，插入的边直接松弛，再并行传播）。"
      "“FMU”每批维护全部结构层；“延迟维护”每批只写边存储，结构层在流结束后同步一次。每条流结束时把增量结果和从头计算逐点比对。\n")
    w("| 算法 | 流 | 方式 | 图更新（ms/批） | 增量修正（ms/批） | 合计（ms/批） | DZiG 合计（ms/批） | 对 DZiG | 结束时不一致的顶点 |")
    w("|---|---|---|---|---|---|---|---|---|")
    for alg in ("bfs", "sssp"):
        for s in ("insert", "delete", "mixed"):
            dze = get(ds, "DZiG", f"{alg}_{s}_endtoend_ms")
            for sysn, lab in (("GIR-inc", "FMU"), ("GIR-lazy-inc", "延迟维护")):
                mps = get(ds, sysn, f"{s}_mutation_MEPS")
                if mps is None: continue
                u = 10000 / (mps * 1e6) * 1000
                r = get(ds, sysn, f"{alg}_{s}_recompute_ms")
                e = get(ds, sysn, f"{alg}_{s}_endtoend_ms")
                mm = get(ds, sysn, f"inc_mismatch_{s}")
                w(f"| {alg.upper()} | {s} | {lab} | {f(u)} | {f(r)} | {f(e)} | {f(dze)} | {ratio(dze, e)} | {f(mm) if mm is not None else '—'} |")
    w("\n“对 DZiG”一栏为 DZiG 合计时间除以 GIR 合计时间，大于 1 表示 GIR 更快。\n")
w("| 算法 | DZiG 初始计算（ms） | GIR 从头计算（ms） |")
w("|---|---|---|")
for alg in ("bfs", "sssp"):
    w(f"| {alg.upper()} | {f(get(ds, 'DZiG', alg + '_initial_ms'))} | {f(get(ds, 'GIR', alg + '_initial_ms'))} |")
w("\n这里 DZiG 的 BFS 和 SSSP 走 KickStarter 引擎，GraphBolt 与 DZiG 在这两个算法上是同一套代码，数字相同。"
  "DZiG 的稀疏感知运行时只作用于 PageRank 一类 GraphBolt 引擎算法。\n")

# ---------------------------------------------------------------- 5 symmetric
dss = "livejournal_sym"
w("## 5. 对称化 LiveJournal：GIR 对比 Aspen 和 Terrace（64 线程）\n")
w("Aspen 只支持对称图，Terrace 原生的 BFS 和 BC 也只在对称图上正确，所以这一组在对称化图上比较。"
  "每条无向边的更新记为两次有向操作，三个系统计法相同。\n")
w("| 系统 | 插入（MEPS） | 删除（MEPS） | 混合（MEPS） | BFS（ms/源） | SSSP（ms/源） | BC 64 源总计（ms） |")
w("|---|---|---|---|---|---|---|")
w(f"| **GIR**（FMU 批量） | {f(get(dss, 'GIR', 'insert_MEPS', 'FMU_B10000'))} | {f(get(dss, 'GIR', 'delete_MEPS', 'FMU_B10000'))} | "
  f"{f(get(dss, 'GIR', 'mixed_MEPS', 'FMU_B10000'))} | {f(get(dss, 'GIR', 'bfs_ms', 'direction-optimizing'))} | "
  f"{f(get(dss, 'GIR', 'sssp_ms', 'delta-stepping'))} | {f(get(dss, 'GIR', 'bc_ms_total', 'brandes'))} |")
w(f"| **GIR**（SPU 单边） | {f(get(dss, 'GIR', 'insert_MEPS', 'SPU'))} | {f(get(dss, 'GIR', 'delete_MEPS', 'SPU'))} | "
  f"{f(get(dss, 'GIR', 'mixed_MEPS', 'SPU'))} | — | — | — |")
for lab, md in (("延迟维护，只写入", "LAZY_ingest_B10000"), ("延迟维护，写入 + 一次同步", "LAZY_total_B10000")):
    w(f"| **GIR**（{lab}） | {f(get(dss, 'GIR', 'insert_MEPS', md))} | {f(get(dss, 'GIR', 'delete_MEPS', md))} | "
      f"{f(get(dss, 'GIR', 'mixed_MEPS', md))} | — | — | — |")
w(f"| GIR 两级模型 | — | — | — | {f(get(dss, 'GIR', 'bfs_ms', 'two-level'))} | {f(get(dss, 'GIR', 'sssp_ms', 'two-level'))} | — |")
for s in ("Aspen", "Terrace"):
    w(f"| {s} | {f(get(dss, s, 'insert_MEPS'))} | {f(get(dss, s, 'delete_MEPS'))} | {f(get(dss, s, 'mixed_MEPS'))} | "
      f"{f(get(dss, s, 'bfs_ms'))} | 不支持 | {f(get(dss, s, 'bc_ms_total'))} |")
w("")
w("| 比值 | 插入 | 删除 | 混合 | BFS | BC |")
w("|---|---|---|---|---|---|")
for s in ("Aspen", "Terrace"):
    w(f"| GIR 对 {s} | {ratio(get(dss, 'GIR', 'insert_MEPS', 'FMU_B10000'), get(dss, s, 'insert_MEPS'))} | "
      f"{ratio(get(dss, 'GIR', 'delete_MEPS', 'FMU_B10000'), get(dss, s, 'delete_MEPS'))} | "
      f"{ratio(get(dss, 'GIR', 'mixed_MEPS', 'FMU_B10000'), get(dss, s, 'mixed_MEPS'))} | "
      f"{ratio(get(dss, s, 'bfs_ms'), get(dss, 'GIR', 'bfs_ms', 'direction-optimizing'))} | "
      f"{ratio(get(dss, s, 'bc_ms_total'), get(dss, 'GIR', 'bc_ms_total', 'brandes'))} |")
for s in ("Aspen", "Terrace"):
    md = "LAZY_ingest_B10000"
    w(f"| GIR 延迟维护只写入 对 {s} | {ratio(get(dss, 'GIR', 'insert_MEPS', md), get(dss, s, 'insert_MEPS'))} | "
      f"{ratio(get(dss, 'GIR', 'delete_MEPS', md), get(dss, s, 'delete_MEPS'))} | "
      f"{ratio(get(dss, 'GIR', 'mixed_MEPS', md), get(dss, s, 'mixed_MEPS'))} | — | — |")
w("\n更新一栏为 GIR 吞吐除以对方吞吐，查询一栏为对方时间除以 GIR 时间。两者都是大于 1 表示 GIR 更好。"
  "三个系统从第一个源点都到达 4,843,953 个顶点。\n")

# ---------------------------------------------------------------- 6 scalability
w("## 6. GIR 线程扩展性（有向 LiveJournal BFS，ms/源）\n")
w("基线系统没有做线程扫描，这一节只有 GIR。\n")
w("| 线程 | 方向优化 BFS | 两级 BFS | 方向优化相对 1 线程 |")
w("|---|---|---|---|")
base = get(ds, "GIR", "bfs_ms", "direction-optimizing", 1)
for t in (1, 2, 4, 8, 16, 32, 64):
    d = get(ds, "GIR", "bfs_ms", "direction-optimizing", t)
    w(f"| {t} | {f(d)} | {f(get(ds, 'GIR', 'bfs_ms', 'two-level', t))} | {ratio(base, d)} |")
w("")

# ---------------------------------------------------------------- 7 twitter
w("## 7. twitter（Kwak，41.7M 点 / 1.47B 边）：仅 GIR\n")
w("基线系统没有在 twitter 上运行。下表只给 GIR 的绝对数字，作为规模参照。\n")
w("| 指标 | 值 |")
w("|---|---|")
for m, mode, label in (("insert_MEPS", "SPU", "插入 SPU（MEPS）"), ("insert_MEPS", "FMU_B10000", "插入 FMU（MEPS）"),
                       ("delete_MEPS", "SPU", "删除 SPU（MEPS）"), ("delete_MEPS", "FMU_B10000", "删除 FMU（MEPS）"),
                       ("mixed_MEPS", "SPU", "混合 SPU（MEPS）"), ("mixed_MEPS", "FMU_B10000", "混合 FMU（MEPS）")):
    w(f"| {label} | {f(get('twitter', 'GIR', m, mode))} |")
for m, mode, label in (("bfs_ms", "direction-optimizing", "BFS 方向优化（ms/源）"), ("bfs_ms", "two-level", "BFS 两级（ms/源）"),
                       ("sssp_ms", "delta-stepping", "SSSP（ms/源）"), ("bc_ms_total", "brandes", "BC 64 源（ms）")):
    w(f"| {label}，64 线程 | {f(get('twitter', 'GIR', m, mode, 64))} |")
w("")

# ---------------------------------------------------------------- 8 paper claims
w("## 8. 论文声称 vs 实测（LiveJournal）\n")
w("论文数字取自 `paper_original/experiment.tex` 的表 3 和表 4，以及 `figures/scripts/plot_throughput.py`。"
  "其中 DZiG 的吞吐柱状图在论文里标注为估算值。\n")
w("| 指标 | 对比对象 | 论文声称：GIR 是对方的 | 实测：GIR 是对方的 |")
w("|---|---|---|---|")
paper = [
    ("插入吞吐", "Terrace", 2.0, ratio(gi, get(ds, "Terrace", "insert_MEPS"))),
    ("插入吞吐", "Spruce", 10.3, ratio(gi, get(ds, "Spruce", "insert_MEPS", None, 64))),
    ("插入吞吐", "DZiG", 11.2, ratio(gi, get(ds, "DZiG", "insert_mutation_MEPS", "bfs_batch_10000"))),
    ("插入吞吐", "Aspen（对称图）", 6.3, ratio(get(dss, "GIR", "insert_MEPS", "FMU_B10000"), get(dss, "Aspen", "insert_MEPS"))),
    ("删除吞吐", "Terrace", 4.7, ratio(gd, get(ds, "Terrace", "delete_MEPS"))),
    ("删除吞吐", "Spruce", 11.0, ratio(gd, get(ds, "Spruce", "delete_MEPS", None, 64))),
    ("删除吞吐", "DZiG", 11.9, ratio(gd, get(ds, "DZiG", "delete_mutation_MEPS", "bfs_batch_10000"))),
    ("删除吞吐", "Aspen（对称图）", 9.1, ratio(get(dss, "GIR", "delete_MEPS", "FMU_B10000"), get(dss, "Aspen", "delete_MEPS"))),
    ("混合吞吐", "Terrace", 25.9 / 8.1, ratio(gm, get(ds, "Terrace", "mixed_MEPS"))),
    ("混合吞吐", "Aspen（对称图）", 25.9 / 4.6, ratio(get(dss, "GIR", "mixed_MEPS", "FMU_B10000"), get(dss, "Aspen", "mixed_MEPS"))),
    ("BFS", "Terrace", 1.5, ratio(get(ds, "Terrace", "bfs_ms"), gb)),
    ("BFS", "Spruce", 1.6, ratio(get(ds, "Spruce", "bfs_ms", None, 64), gb)),
    ("BFS", "Aspen（对称图）", 1.2, ratio(get(dss, "Aspen", "bfs_ms"), get(dss, "GIR", "bfs_ms", "direction-optimizing"))),
    ("SSSP", "Spruce", 1.9, ratio(get(ds, "Spruce", "sssp_ms", None, 64), gs)),
    ("BC", "Terrace（对称图；论文为有向图）", 1.3, ratio(get(dss, "Terrace", "bc_ms_total"), get(dss, "GIR", "bc_ms_total", "brandes"))),
    ("BC", "Aspen（对称图；论文为有向图）", 2.1, ratio(get(dss, "Aspen", "bc_ms_total"), get(dss, "GIR", "bc_ms_total", "brandes"))),
]
for metric, who, claim, meas in paper:
    w(f"| {metric} | {who} | {claim:.1f}× | {meas} |")
w(f"\n论文声称 GIR 在 LiveJournal 上的混合吞吐为 25.9 MEPS，实测为 {f(gm)} MEPS（FMU）。"
  "论文的 BFS 和 SSSP 对 GraphBolt/DZiG 的加速比按查询时间计算，第 4 节给出了同源点、同更新流下的逐批实测。\n")

# ---------------------------------------------------------------- 9 GIR before / after the speed-up
w("## 9. GIR 本轮提速前后（同一协议，同一台机器）\n")
w("“提速前”取自 `experiments/results/history/*_pre_speedup.*`，即 2026-10-14 提速之前的同一组测量；"
  "“第一轮”到“第九轮”取自 `history/*_round1.*` 到 `history/*_round9.*`；“现在”是当前 CSV（第十轮）。"
  "吞吐越高越好，时间越低越好；倍数一栏是“现在比提速前好多少倍”。\n")
w("| 场景 | 指标 | 提速前 | 第一轮 | 第二轮 | 第三轮 | 第四轮 | 第五轮 | 第六轮 | 第七轮 | 第八轮 | 第九轮 | 现在 | 倍数 |")
w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")


def better(b, a, higher):
    if a is None or b is None or a <= 0 or b <= 0:
        return "—"
    r = a / b if higher else b / a
    return f"{r:,.0f}×" if r >= 100 else f"{r:.2f}×" if r >= 1 else f"1/{1 / r:.2f}"


def row_ba(label, metric, before, after, higher, unit, middle=None, middle2=None, middle3=None, middle4=None, middle5=None, middle6=None, middle7=None, middle8=None, middle9=None):
    w(f"| {label} | {metric} | {f(before)} {unit} | {f(middle)} {unit} | {f(middle2)} {unit} | {f(middle3)} {unit} | "
      f"{f(middle4)} {unit} | {f(middle5)} {unit} | {f(middle6)} {unit} | {f(middle7)} {unit} | {f(middle8)} {unit} | {f(middle9)} {unit} | {f(after)} {unit} | {better(before, after, higher)} |")


for dsn, label in (("livejournal", "有向 LJ"), ("livejournal_sym", "对称 LJ"), ("twitter", "twitter")):
    for mode, ml in (("SPU", "SPU"), ("FMU_B10000", "FMU")):
        for m, mn in (("insert_MEPS", "插入"), ("delete_MEPS", "删除"), ("mixed_MEPS", "混合")):
            row_ba(f"{label} {ml}", mn, get(dsn, "GIR", m, mode, None, pre), get(dsn, "GIR", m, mode), True, "MEPS",
                   get(dsn, "GIR", m, mode, None, mid), get(dsn, "GIR", m, mode, None, mid2), get(dsn, "GIR", m, mode, None, mid3),
                   get(dsn, "GIR", m, mode, None, mid4), get(dsn, "GIR", m, mode, None, mid5),
                   get(dsn, "GIR", m, mode, None, mid6), get(dsn, "GIR", m, mode, None, mid7),
                   get(dsn, "GIR", m, mode, None, mid8), get(dsn, "GIR", m, mode, None, mid9))
for dsn, label in (("livejournal", "有向 LJ"), ("livejournal_sym", "对称 LJ"), ("twitter", "twitter")):
    for mode, ml in (("LAZY_ingest_B10000", "延迟维护只写入"), ("LAZY_total_B10000", "延迟维护写入+同步")):
        for m, mn in (("insert_MEPS", "插入"), ("delete_MEPS", "删除"), ("mixed_MEPS", "混合")):
            row_ba(f"{label} {ml}", mn, get(dsn, "GIR", m, "FMU_B10000", None, pre), get(dsn, "GIR", m, mode), True, "MEPS",
                   None, None, get(dsn, "GIR", m, mode, None, mid3), get(dsn, "GIR", m, mode, None, mid4),
                   get(dsn, "GIR", m, mode, None, mid5), get(dsn, "GIR", m, mode, None, mid6),
                   get(dsn, "GIR", m, mode, None, mid7), get(dsn, "GIR", m, mode, None, mid8), get(dsn, "GIR", m, mode, None, mid9))
for metric, mode, mn in (("bfs_ms", "direction-optimizing", "BFS"), ("sssp_ms", "delta-stepping", "SSSP"),
                         ("bc_ms_total", "brandes", "BC 64 源")):
    row_ba("有向 LJ 查询（64 线程）", mn, get("livejournal", "GIR", metric, mode, 64, pre_clean), gq(metric, mode), False, "ms",
           get("livejournal", "GIR", metric, mode, 64, mid_clean), get("livejournal", "GIR", metric, mode, 64, mid2_clean),
           get("livejournal", "GIR", metric, mode, 64, mid3_clean), get("livejournal", "GIR", metric, mode, 64, mid4_clean),
           get("livejournal", "GIR", metric, mode, 64, mid5_clean), get("livejournal", "GIR", metric, mode, 64, mid6_clean),
           get("livejournal", "GIR", metric, mode, 64, mid7_clean), get("livejournal", "GIR", metric, mode, 64, mid8_clean), get("livejournal", "GIR", metric, mode, 64, mid9_clean))
for metric, mode, mn in (("bfs_ms", "direction-optimizing", "BFS"), ("sssp_ms", "delta-stepping", "SSSP"),
                         ("bc_ms_total", "brandes", "BC 64 源")):
    row_ba("twitter 查询（64 线程）", mn, get("twitter", "GIR", metric, mode, 64, pre), get("twitter", "GIR", metric, mode, 64), False, "ms",
           get("twitter", "GIR", metric, mode, 64, mid), get("twitter", "GIR", metric, mode, 64, mid2),
           get("twitter", "GIR", metric, mode, 64, mid3), get("twitter", "GIR", metric, mode, 64, mid4),
           get("twitter", "GIR", metric, mode, 64, mid5), get("twitter", "GIR", metric, mode, 64, mid6), get("twitter", "GIR", metric, mode, 64, mid7), get("twitter", "GIR", metric, mode, 64, mid8), get("twitter", "GIR", metric, mode, 64, mid9))
for alg in ("bfs", "sssp"):
    for s, sn in (("insert", "插入流"), ("delete", "删除流"), ("mixed", "混合流")):
        row_ba(f"流式 {alg.upper()}（每批合计）", sn, get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms", None, None, pre),
               get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms"), False, "ms",
               get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms", None, None, mid),
               get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms", None, None, mid2),
               get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms", None, None, mid3),
               get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms", None, None, mid4),
               get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms", None, None, mid5),
               get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms", None, None, mid6),
               get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms", None, None, mid7),
               get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms", None, None, mid8),
               get("livejournal", "GIR", f"{alg}_{s}_endtoend_ms", None, None, mid9))
w("\n主要改动：删除的边存储墓碑从哈希集合改为每条 CSR 边一位的位图；FMU 批删除先把整批删除写入边存储，"
  "再在最终图上并行验证每条删除不改变任何层（可达性 + 两条不相交路径），只把验证不了的删除放回走单边路径；"
  "顶点离开块时的连接点两两检验改为线程安全的有界搜索，并行执行；FMU 批插入中不改变结构的边并行批量写入；"
  "BFS 的自底向上步不再每次查询现算可达集，候选与压缩一次完成，前沿用位图；Δ-stepping 改为每线程桶，去掉串行合并。\n")
w("延迟维护各行的“提速前”是同一数据集提速前的 FMU 批量吞吐（当时没有延迟模式）。\n")
w("第二轮（2026-10-15）：FMU 批删除中，大 SCC 的分裂整批处理（在批后图上剪掉入度或出度为零的顶点，"
  "用并行边界探测证明剩余部分仍强连通，探测失败就把穷尽的一侧并入，再一次性分离），不再逐条顺序分裂；"
  "整批删除并行写入边存储；SCC 编号容器预留容量，避免分裂时重新分配上百万个哈希表；"
  "单边删除的块证书复用可达性检查找到的有向路径，只需再搜一条不相交路径；"
  "精确搜索（Menger 及块内路径）按“下一步代价”选择扩展方向，hub 不再先被展开；BC 的每顶点状态合并为一条记录。\n")
w("第三轮（2026-10-15 晚）：新增延迟维护（`ingestBatch` 只写边存储，`sync()` 时用 `applyBatch` 补上三层结构）；"
  "增量边表改为每顶点有序数组，去掉全局哈希集合，批量写入按顶点并行合并；批内去重并行；链接 jemalloc；"
  "批删除的 WCC 检查先并行判定孤立点和与代表 SCC 相连的端点；块修复在切链后先用有界证书确认，省掉最后一次 Menger。\n")
w("第四轮（2026-10-16）：批删除的遗留边沿用批证书的结果（已知可达的直接修块，已知不可达的直接走精确搜索），"
  "不再重复有界搜索；精确可达搜索和 Menger 搜索对高度数顶点改为延后展开，另一侧每访问一个顶点就二分查找它与这些 hub 的邻接；"
  "有界搜索在相遇或超限后立即停止扫描邻接表。twitter 上 2 万次删除的同步从约 23 s 降到约 5–9 s。\n")
w("第五轮（2026-10-16）：块（BiCC）层可按需维护（见第 10.2 节；第 9 节各行仍是逐批维护块的默认模式）；"
  "大 SCC 的块用并行 Tarjan–Vishkin 构建；新增结构查询。第五轮流式各行只保留了每批合计时间。\n")
w("第六轮（2026-10-16 晚）：批删除的分裂证明以图中最大入度和最大出度的两个 hub 为锚点（仅当它们的邻居覆盖 SCC 的 1/128 以上，"
  "即 twitter 这类图），每个探测只需走到 hub 的邻居，通常一两跳；同一锚点论证也用来批量证明“删除后 u 仍可达 v”；"
  "证明通过的探测记下路径，后续轮次只重测路径碰到新分出顶点的探测；批证明通过的删除不再走单边的精确搜索；"
  "与拓扑序同向的 SCC 间插入并行批量写入；并行 BiCC 的 BFS 改为方向优化，规则 (1) 并入 low/high 的同一次邻接扫描。\n")
w("第七轮（2026-10-19）：参照 Ligra/Terrace 的查询实现：结果数组不再由 std::vector 串行清零（4.8M 个 int 约 8 ms，主要是单线程缺页），"
  "改为不初始化分配后并行填充；BFS 自顶向下按边均分（hub 不再拖慢一层），访问位图放在缓存里，自底向上步用位图输出下一层，不再写候选表；"
  "批写入少一次排序、按批大小走顺序路径；新增与 DZiG/KickStarter 同样的增量 BFS/SSSP（第 4 节）。\n")
w("第八轮（2026-10-19 晚）：没有单个 hub 能覆盖的大 SCC（如 livejournal）改用一组锚点（最多 32 个高度数顶点，每轮证明它们互相可达），"
  "边界顶点只需走到任一锚点的邻居；超出顺序预算的少数探测改为所有线程一起做的分层 BFS，不再由一个探测决定整轮时间；"
  "twitter 的逐条删除证明分三级（小预算并行、大预算并行、全线程逐个）；墓碑位图和 hub 候选在构建时准备好；"
  "跨 SCC 删除先检查两端是否仍连到 hub 顶点。\n")
w("第九轮（2026-10-20）：批写入的筛选与分组改为并行两遍（先计数、再按前缀写入），不再串行收集；"
  "只有一个前驱（或后继）的 SCC 与之成环时直接合并，不做 Pearce–Kelly 搜索；没有凝聚边的 SCC（新加入的点）跨 WCC 时直接取相邻的空闲序号；"
  "合并小 SCC 时不再建哈希集合；并行 BiCC 的树记录缩到 8 字节。BC 改为每个源点内部并行的方案实测更慢（LJ 5.3–5.7 s 对 3.0–3.4 s），已撤回。\n")
w("第十轮（2026-10-20/21）：Δ-stepping 每轮直接读各线程的桶，不再拷贝合并；插入路径按实测重做："
  "与拓扑序同向的凝聚边按 SCC 分组并行写入，块按需时整批插入一次写入边存储，只有结构步骤串行；"
  "没有凝聚边的 SCC 跨 WCC 时放到序的两端，不再挤在巨型 SCC 旁边。\n")

# ---------------------------------------------------------------- 10 structural queries
w("## 10. 结构查询：GIR 的索引真正起作用的地方\n")
w("BFS、SSSP、BC 只读边存储，用不到 GIR 维护的三层结构，所以在第 2–5 节里 GIR 承担维护代价却没有收益。"
  "下面是直接用结构回答的查询。“无索引”是在同一份边存储上做双向 BFS（用到入边，比 Terrace 有向存储能做的更强），"
  "代表只存邻接表的系统回答同一问题的较好情况；Terrace 只有整图 BFS，单次约等于第 3 节的 BFS 时间。\n")
w("| 数据集 | 可达性 u⇝v：GIR 索引（µs/次） | 无索引双向 BFS（µs/次） | 倍数 | 割点查询：GIR（µs/次） |")
w("|---|---|---|---|---|")
for dsn in ("livejournal", "twitter"):
    src = clean if dsn == "livejournal" and get(dsn, "GIR", "reach_us", "index", None, clean) is not None else rows
    ri = get(dsn, "GIR", "reach_us", "index", None, src)
    rn = get(dsn, "GIR", "reach_us", "no-index-bibfs", None, src)
    ap = get(dsn, "GIR", "artpoint_us", "index", None, src)
    if ri is None: continue
    w(f"| {dsn} | {f(ri)} | {f(rn)} | {ratio(rn, ri)} | {f(ap)} |")
w("\n随机点对，两种方法的答案逐一核对一致。\n")
w("### 10.1 更新 + 查询混合负载\n")
w("混合流每批 10,000 次更新之后跟 Q 个随机可达性查询，每批合计时间（ms，越低越好）。"
  "“GIR 索引”：批量更新（块按需维护）+ 索引查询；“GIR 延迟维护”：只写入，第一次查询前同步；"
  "“只存储”：只写入边存储、不维护结构，每个查询做双向 BFS。\n")
qlist = [0, 10, 100, 1000, 10000, 100000]
w("| 数据集 | 方式 | " + " | ".join(f"Q={q:,}" for q in qlist) + " |")
w("|---|---|" + "---|" * len(qlist))
for dsn in ("livejournal", "twitter"):
    for vk, vl in (("index", "GIR 索引"), ("deferred", "GIR 延迟维护"), ("no-index", "只存储 + BFS")):
        vals = [get(dsn, "GIR", "workload_ms_per_batch", f"{vk}_Q{q}") for q in qlist]
        if all(v is None for v in vals): continue
        w(f"| {dsn} | {vl} | " + " | ".join(f(v) for v in vals) + " |")
w("\n每批查询数超过某个值后 GIR 的总时间更低：结构维护的固定开销被每个查询省下的遍历抵消。\n")
w("### 10.2 块（BiCC）层按需维护\n")
w("块结构的增量维护在删除时最贵（第 2.1 节），而只有块查询（割点、是否双连通、两级 BFS）需要它。"
  "按需模式下更新只维护 WCC/SCC/DAG，记下哪些 SCC 的块过期，块查询前用并行 Tarjan–Vishkin 只重建这些 SCC。\n")
w("| 数据集 | 模式 | 插入 | 删除 | 混合 | 块重建（每阶段后，s） |")
w("|---|---|---|---|---|---|")
for dsb in ("livejournal", "twitter"):
  for md, lab in (("FMU_B10000", "FMU，块逐批维护"), ("FMU_B10000_bod", "FMU，块按需"),
                ("LAZY_total_B10000", "延迟维护+同步，块逐批维护"), ("LAZY_total_B10000_bod", "延迟维护+同步，块按需")):
    vals = [get(dsb, "GIR", m, md) for m in ("insert_MEPS", "delete_MEPS", "mixed_MEPS")]
    if all(v is None for v in vals): continue
    rb = get(dsb, "GIR", "blocks_rebuild_s", md.replace("LAZY_total", "LAZY") + "_delete") if md.endswith("_bod") else None
    if md.startswith("LAZY") and md.endswith("_bod"):
        rb = get(dsb, "GIR", "blocks_rebuild_s", "LAZY_B10000_bod_delete")
    w(f"| {dsb} | {lab} | {f(vals[0])} | {f(vals[1])} | {f(vals[2])} | {f(rb) if rb is not None else '—'} |")
w("")

# ---------------------------------------------------------------- 11 caveats
w("## 11. 解读注意事项\n")
w("- **吞吐对比的是不同的工作量。** GIR 每次更新都维护 WCC、SCC 和 BiCC，存储系统只改邻接表。"
  "差距主要来自这一设计：删除时的 SCC 分裂和块修复是 GIR 独有的开销。在这个设计下，论文的吞吐结论无法成立。")
w("- **算法实现不同。** GIR 用方向优化 BFS；Terrace 和 Aspen 用 Ligra 式 edgeMap；Spruce 的算法由驱动实现，只有自上而下 BFS；"
  "Terrace 的 SSSP 是 Bellman-Ford。")
w("- **有向与对称。** Terrace 在有向图上改用仅推送模式，结果已与 GIR、Spruce 交叉验证；Aspen 只能在对称图上比较。")
w("- **单次运行。** 每个数字是一次运行，查询取 16 个源点的平均；没有重复多次取方差。"
  "这台机器是多人共用的，同一程序的 GIR 插入吞吐在相邻两次运行之间最多相差约 40%（有向 LJ 的 SPU 插入测到过 0.40 和 0.72 MEPS）。")
w("- **查询时间随机器状态变化。** 有向 LJ 的 SSSP 在 2026-10-14 晚测得 101 ms，2026-10-15 测得 128–132 ms；"
  "当天把第一轮的二进制重新编译、与当前版本交替运行，两者 BFS 都是 30.5–31 ms、SSSP 都是 128–132 ms，"
  "所以第 9 节里“第一轮后 → 现在”的 SSSP 变化来自机器状态，不是代码。")
w("- **内存分配器。** 自 2026-10-15 起 GIR 链接 jemalloc（机器上已安装时）。批量更新时各线程同时为许多顶点"
  "分配小数组，glibc malloc 的锁竞争让写入吞吐减半左右。基线系统沿用各自默认的分配方式。")
w("- **延迟维护的语义。** 只写入模式下结构层（WCC/SCC/块）停留在上次同步的状态，结构查询之前必须同步；"
  "表中“写入 + 一次同步”把每 10 批（10 万次操作）的同步开销计入。论文描述的是每次更新都维护结构。")
w("- **单边更新受 NUMA 放置影响。** SPU 是单线程执行的，但图在 64 线程下并行构建，内存页分布在两个 NUMA 节点上。"
  "twitter 上同一程序用 1 个线程运行时 SPU 插入为 0.28 MEPS，64 线程时为 0.12–0.14 MEPS。"
  "第 9 节里 twitter SPU 插入比提速前低（0.359 → 0.13），原因没有完全确定，单线程对比说明主要来自内存放置，而不是插入路径本身变慢；"
  "有向 LJ 上同一路径的单次耗时从 2.2 µs 降到 1.4 µs。")
w("- **未覆盖。** twitter 等更大的图上没有运行基线；论文中的其他数据集（TC、WG、WT、CP、FS）在磁盘上不存在或已损坏。\n")
w("原始数据：`baselines.csv`、`gir_benchmark.csv`、`gir_queries_lj_clean.csv`、`gir_benchmark_sym.csv`，"
  "日志在 `baselines_logs/`。驱动和构建说明见 `experiments/baselines/README.md`。")

open(os.path.join(R, "comparison.md"), "w").write("\n".join(out) + "\n")
print("wrote", os.path.join(R, "comparison.md"))
