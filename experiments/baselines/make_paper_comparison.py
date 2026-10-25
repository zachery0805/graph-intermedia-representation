#!/usr/bin/env python3
"""Paper-style comparison (Sec. 7 of the paper) from experiments/results/paper/<dataset>/.

Every experiment of the paper's evaluation, all seven datasets, every compared
system, next to the numbers the paper claims. Output:
experiments/results/paper_comparison.md.
"""
import csv, os, re, math

R = "/home/cxl6029/projects/GIR/experiments/results"
PD = os.environ.get("PAPER_DIR", os.path.join(R, "paper"))
OUTF = os.environ.get("PAPER_OUT", os.path.join(R, "paper_comparison.md"))
DATASETS = [("twitter_combined", "TC"), ("web-Google", "WG"), ("wiki-topcats", "WT"), ("cit-Patents", "CP"),
            ("livejournal", "LJ"), ("twitter", "TW"), ("friendster", "FS")]
ABBR = [a for _, a in DATASETS]

# ---------------------------------------------------------------- paper claims
P_INS = {"Terrace": [1.7, 1.5, 4.0, 3.9, 2.0, 2.8, 2.7, 2.7], "Aspen": [1.3, 1.2, 6.5, 2.8, 6.3, 3.8, 2.0, 3.4],
         "DZiG": [9.4, 9.6, 9.6, 17.2, 11.2, 6.7, 5.6, 9.9], "Spruce": [8.7, 8.7, 9.0, 16.2, 10.3, 6.0, 5.1, 9.1]}
P_DEL = {"Terrace": [2.9, 2.1, 4.6, 3.3, 4.7, 4.0, 2.5, 3.4], "Aspen": [2.5, 1.7, 7.8, 3.8, 9.1, 3.5, 2.2, 4.4],
         "DZiG": [8.0, 9.1, 10.2, 14.7, 11.9, 6.5, 4.6, 9.3], "Spruce": [7.3, 8.4, 9.5, 13.5, 11.0, 5.8, 4.0, 8.5]}
P_MIX = {"Spruce": [3.5, 3, 3, 3.1, 2.7, 2.9, 2.5], "Aspen": [20.2, 18.9, 4.3, 31.2, 4.6, 8.9, 10],
         "Terrace": [17.5, 16.4, 6.7, 31.2, 8.1, 7.6, 7.8], "GIR": [26.3, 24.6, 27.6, 43.7, 25.9, 16, 14]}
P_ALG = {("BFS", "Terrace"): [1.5, 1.2, 1.2, 1.8, 1.5, 1.1, 1.1, 1.3], ("BFS", "Aspen"): [1.1, 1.2, 1.4, 1.4, 1.2, 1.2, 1.6, 1.3],
         ("BFS", "Spruce"): [3, 2.7, 1.4, 8.6, 1.6, 1.8, 1.9, 3], ("BFS", "GraphBolt"): [1.4, 1.3, 1.2, 1.7, 1.3, 1.1, 1.1, 1.3],
         ("BFS", "DZiG"): [1.3, 1.2, 1.1, 1.5, 1.2, 1.1, 1.1, 1.2],
         ("SSSP", "Spruce"): [7, 2.4, 3.1, 1.9, 1.9, 1.1, 1.4, 2.7], ("SSSP", "GraphBolt"): [2.0, 1.8, 1.8, 2.2, 1.5, 1.2, 1, 1.2],
         ("SSSP", "DZiG"): [1.7, 1.6, 1.7, 1.9, 1.3, 1.1, 1.2, 1.6],
         ("BC", "Terrace"): [2.1, 1.8, 1.9, 1.5, 1.3, 1.5, 1.5, 1.7], ("BC", "Aspen"): [2, 2.1, 2.1, 1.6, 2.1, 1.2, 1.1, 1.7]}
P_SCALE = {"GIR": {"TC": 11.9, "WG": 7.5, "WT": 47.4, "CP": 11, "LJ": 45, "TW": 8.9, "FS": 11.1, "Avg": 20.9},
           "Spruce": {"WT": 17, "Avg": 17}, "Aspen": {"WT": 43.6}}
P_DS = {"TC": ("81.3K", "2.4M"), "WG": ("876K", "5.1M"), "WT": ("1.8M", "28.5M"), "CP": ("3.8M", "16.5M"),
        "LJ": ("4.9M", "69M"), "TW": ("61.6M", "2.4B"), "FS": ("65.6M", "1.8B")}

# ---------------------------------------------------------------- data
rows = []  # (dataset, system, metric, mode, value, threads)


def fnum(x):
    try:
        return float(x)
    except ValueError:
        return None


for ds, _ in DATASETS:
    d = os.path.join(PD, ds)
    if not os.path.isdir(d):
        continue
    for fn, sysname in (("gir.csv", "GIR"), ("gir_bod.csv", "GIR-bod"), ("gir_sym.csv", "GIR")):
        p = os.path.join(d, fn)
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p)):
            v = fnum(r["value"])
            if v is None:
                continue
            mode = r["mode"]
            s = sysname
            if fn == "gir_bod.csv":
                s = "GIR"
            rows.append((r["dataset"], s, r["metric"], mode, v, int(r["threads"])))
    p = os.path.join(d, "baselines.csv")
    parts = [p] if os.path.exists(p) else sorted(os.path.join(d, x) for x in os.listdir(d) if x.startswith("b_") and x.endswith(".csv"))
    for p in parts:
        for line in open(p):
            f = line.strip().split(",")
            if len(f) < 8 or f[0] == "dataset":
                continue
            v = fnum(f[4])
            rows.append((f[0], f[1], f[2], f[3], v if v is not None else float("nan"), int(f[7]) if f[7].isdigit() else 0))


def get(ds, system, metric, mode=None, threads=None, pick=max):
    """Value of the matching row; with several (threads not given) the one at the most threads."""
    c = [(t, v) for (d, s, m, mo, v, t) in rows if d == ds and s == system and m == metric
         and (mode is None or mo == mode or (mode.endswith("*") and mo.startswith(mode[:-1])))
         and (threads is None or t == threads)]
    if not c:
        # a recorded failure (mode "FAILED:...") shows as "失败"
        if any(d == ds and s == system and m == metric and mo.startswith("FAILED") for (d, s, m, mo, v, t) in rows):
            return -1.0
        return None
    t = pick(x[0] for x in c)
    return [v for (tt, v) in c if tt == t][-1]


def ths(ds, system, metric, mode):
    return sorted({t for (d, s, m, mo, v, t) in rows if d == ds and s == system and m == metric and mo == mode})


def f(x, unit=""):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    if x < 0:
        return "失败"
    a = abs(x)
    s = f"{x:,.0f}" if a >= 1000 else f"{x:.0f}" if a >= 100 else f"{x:.1f}" if a >= 10 else f"{x:.2f}" if a >= 1 else f"{x:.3g}"
    return s + unit


def ratio(a, b):
    """a / b, as '3.2×' or '1/3.2'."""
    if a is None or b is None or a <= 0 or b <= 0 or math.isnan(a) or math.isnan(b):
        return "—"
    r = a / b
    return f"{r:.2f}×" if r >= 1 else f"1/{1 / r:.1f}"


def geo(xs):
    xs = [x for x in xs if x is not None and x > 0 and not math.isnan(x)]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else None


def logval(ds, pat, fn="gir.log"):
    p = os.path.join(PD, ds, fn)
    if not os.path.exists(p):
        return None
    m = re.search(pat, open(p).read())
    return m.group(1) if m else None


out = []
w = out.append
w("# GIR 论文实验复现：全部数据集、全部对比系统\n")
w("由 `experiments/baselines/make_paper_comparison.py` 从 `experiments/results/paper/<数据集>/` 生成。"
  "每个数据集用 `experiments/baselines/run_paper.sh`（或 `run_paper2.sh`）依次运行全部系统，同一台机器，任何时刻只跑一个系统。"
  "“论文”列是论文正文、表格或作图脚本里的数字；“实测”列是本次运行。\n")

# ---------------------------------------------------------------- setup
w("## 0. 设置与数据集\n")
w("- 机器：2 × AMD EPYC 9334（64 核 / 128 线程），1.1 TB 内存；默认 64 线程（与论文一致）。")
w("- 协议（论文第 7.2 节）：每个图取一半边建初始图，另一半作为插入流，删除从已加载的一半中抽取；混合流把两者打乱交错。"
  "每个阶段 100,000 次操作、每批 10,000（论文写的是整条留出的一半；这里每阶段取 10 万次，所有系统相同）。")
w("- BFS/SSSP/BC 在完整图上跑，16 个源点取平均（BC 为 64 个源点总时间），所有系统用同一组源点。SSSP 权重为 1–100 的确定性伪随机值。")
w("- 有向与对称：Aspen 只支持对称图，Terrace 的 BC 只在对称图上正确，所以这两项在对称化图上测，并与同一对称图上的 GIR 比较；"
  "Terrace 的 BFS 在有向图上用仅推送版本（其稠密 edgeMap 在有向图上结果不对）。其余系统都在有向图上。")
w("- DZiG 与 GraphBolt 的 BFS/SSSP 是同一个 KickStarter 引擎（同一份代码），一组数字同时代表两者；它们的“从头计算”在 50% 初始图上，"
  "所以与 GIR 在同一 50% 图上的从头计算比较。\n")
w("| 数据集 | 缩写 | 论文 点 / 边 | 实测 点 / 边（去自环、去重） | SCC 数 | WCC 数 | 备注 |")
w("|---|---|---|---|---|---|---|")
for ds, ab in DATASETS:
    V = logval(ds, r"V=(\d+) E=")
    E = logval(ds, r"V=\d+ E=(\d+)")
    scc = logval(ds, r"Found (\d+) SCCs")
    wcc = logval(ds, r"Found (\d+) WCCs")
    note = {"TW": "Kwak 2010 版（41.7M / 1.47B）；论文的 61.6M / 2.4B 与 SNAP 和 Kwak 都对不上",
            "FS": "SNAP 无向表每条边只列一次（u < v），按列出方向是 DAG；这里按端点哈希随机定向（所有系统同一个图）",
            "TC": "论文边数 2.4M 与 SNAP（1.77M）不符"}.get(ab, "")
    vv = f"{int(V):,} / {int(E):,}" if V else "未运行"
    w(f"| {ds} | {ab} | {P_DS[ab][0]} / {P_DS[ab][1]} | {vv} | {int(scc):,} | {int(wcc):,} | {note} |" if scc else
      f"| {ds} | {ab} | {P_DS[ab][0]} / {P_DS[ab][1]} | {vv} | — | — | {note} |")
w("")

# ---------------------------------------------------------------- update throughput
GIR_MODES = [("SPU", "GIR 单边（SPU）"), ("FMU_B10000", "GIR 批量（FMU，块逐批维护）"),
             ("FMU_B10000_bod", "GIR 批量（FMU，块按需）"), ("LAZY_total_B10000_bod", "GIR 延迟维护 + 同步（块按需）"),
             ("LAZY_ingest_B10000", "GIR 只写边存储（延迟维护的写入部分）")]


def upd(ds, system, metric):
    if system == "Terrace":
        return get(ds, "Terrace", metric)
    if system == "Aspen":
        return get(ds + "_sym", "Aspen", metric)
    if system == "Spruce":
        return get(ds, "Spruce", metric)
    if system == "DZiG":
        m = {"insert_MEPS": "insert_mutation_MEPS", "delete_MEPS": "delete_mutation_MEPS", "mixed_MEPS": "mixed_mutation_MEPS"}[metric]
        return get(ds, "DZiG", m, "bfs_batch_10000")
    return None


def throughput_section(title, metric, paper):
    w(f"### {title}\n")
    w("吞吐（MEPS，越高越好）。Aspen 一行是对称图上的吞吐（每条无向更新计两次操作），与之比较的 GIR 也取对称图。\n")
    w("| 系统 | " + " | ".join(ABBR) + " |")
    w("|---|" + "---|" * len(ABBR))
    for mode, lab in GIR_MODES:
        w(f"| {lab} | " + " | ".join(f(get(ds, "GIR", metric, mode)) for ds, _ in DATASETS) + " |")
    w("| GIR 对称图（FMU，块逐批维护） | " + " | ".join(f(get(ds + "_sym", "GIR", metric, "FMU_B10000")) for ds, _ in DATASETS) + " |")
    for sname in ("Terrace", "Aspen", "Spruce", "DZiG"):
        w(f"| {sname} | " + " | ".join(f(upd(ds, sname, metric)) for ds, _ in DATASETS) + " |")
    w("")
    if paper is None:
        return
    w("GIR 相对各系统的倍数（GIR ÷ 对方，大于 1 表示 GIR 更快；1/N 表示 GIR 慢 N 倍）。"
      "“论文”为作图脚本 `paper_original/figures/scripts/plot_throughput.py` 里的数字。实测用 GIR 批量（块按需）；"
      "对 Aspen 用对称图上的 GIR 批量。\n")
    w("| 对比 | 来源 | " + " | ".join(ABBR) + " | 平均 |")
    w("|---|---|" + "---|" * (len(ABBR) + 1))
    for sname in ("Terrace", "Aspen", "DZiG", "Spruce"):
        meas, vals = [], []
        for ds, _ in DATASETS:
            g = get(ds + "_sym", "GIR", metric, "FMU_B10000") if sname == "Aspen" else get(ds, "GIR", metric, "FMU_B10000_bod")
            o = upd(ds, sname, metric)
            meas.append(ratio(g, o))
            vals.append(g / o if g and o and o > 0 and g > 0 else None)
        gm = geo(vals)
        w(f"| GIR ÷ {sname} | 论文 | " + " | ".join(f"{x}×" for x in paper[sname][:7]) + f" | {paper[sname][7]}× |")
        w(f"| GIR ÷ {sname} | 实测 | " + " | ".join(meas) + f" | {ratio(gm, 1.0) if gm else '—'} |")
    w("\n实测平均为几何平均。\n")


w("## 1. 更新吞吐（论文第 7.2 节，图 7、图 8、表 3）\n")
throughput_section("1.1 只插入（论文图 7）", "insert_MEPS", P_INS)
throughput_section("1.2 只删除（论文图 8）", "delete_MEPS", P_DEL)
throughput_section("1.3 混合（论文表 3）", "mixed_MEPS", None)
w("论文表 3 的混合吞吐（MEPS）与实测对照：\n")
w("| 系统 | 来源 | " + " | ".join(ABBR) + " |")
w("|---|---|" + "---|" * len(ABBR))
for sname in ("GIR", "Terrace", "Aspen", "Spruce"):
    w(f"| {sname} | 论文 | " + " | ".join(f(x) for x in P_MIX[sname]) + " |")
    if sname == "GIR":
        w(f"| GIR（FMU，块按需） | 实测 | " + " | ".join(f(get(ds, "GIR", "mixed_MEPS", "FMU_B10000_bod")) for ds, _ in DATASETS) + " |")
        w(f"| GIR（FMU，块逐批维护） | 实测 | " + " | ".join(f(get(ds, "GIR", "mixed_MEPS", "FMU_B10000")) for ds, _ in DATASETS) + " |")
    else:
        w(f"| {sname} | 实测 | " + " | ".join(f(upd(ds, sname, "mixed_MEPS")) for ds, _ in DATASETS) + " |")
w("\nSpruce 的混合流在把某顶点最后一条边删掉后，再往该顶点插边会卡死（已知缺陷，见 comparison.md），超时记为“失败”。\n")

# ---------------------------------------------------------------- analytics
w("## 2. 图分析（论文第 7.3 节，表 4）\n")
w("时间（ms，越低越好）。BFS/SSSP 为每个源点平均，BC 为 64 个源点总计。64 线程。\n")


def alg(ds, system, a):
    if a == "BFS":
        return {"GIR": get(ds, "GIR", "bfs_ms", "direction-optimizing", 64, max) or get(ds, "GIR", "bfs_ms", "direction-optimizing"),
                "GIR 对称图": get(ds + "_sym", "GIR", "bfs_ms", "direction-optimizing"),
                "Terrace": get(ds, "Terrace", "bfs_ms", "push-only", 64) or get(ds, "Terrace", "bfs_ms", "push-only"),
                "Aspen": get(ds + "_sym", "Aspen", "bfs_ms", "stock", 64) or get(ds + "_sym", "Aspen", "bfs_ms", "stock"),
                "Spruce": get(ds, "Spruce", "bfs_ms", "top-down", 64) or get(ds, "Spruce", "bfs_ms", "top-down"),
                "GIR 50% 图": get(ds, "GIR", "bfs_initial_ms", "from-scratch"),
                "DZiG/GraphBolt 50% 图": get(ds, "DZiG", "bfs_initial_ms", "from-scratch")}.get(system)
    if a == "SSSP":
        return {"GIR": get(ds, "GIR", "sssp_ms", "delta-stepping"),
                "Spruce": get(ds, "Spruce", "sssp_ms", "delta-stepping"),
                "GIR 50% 图": get(ds, "GIR", "sssp_initial_ms", "from-scratch"),
                "DZiG/GraphBolt 50% 图": get(ds, "DZiG", "sssp_initial_ms", "from-scratch")}.get(system)
    if a == "BC":
        return {"GIR": get(ds, "GIR", "bc_ms_total", "brandes"),
                "GIR 对称图": get(ds + "_sym", "GIR", "bc_ms_total", "brandes"),
                "Terrace": get(ds + "_sym", "Terrace", "bc_ms_total", "stock"),
                "Aspen": get(ds + "_sym", "Aspen", "bc_ms_total", "stock"),
                "Spruce": get(ds, "Spruce", "bc_ms_total", "brandes")}.get(system)


SYS = {"BFS": ["GIR", "GIR 对称图", "Terrace", "Aspen", "Spruce", "GIR 50% 图", "DZiG/GraphBolt 50% 图"],
       "SSSP": ["GIR", "Spruce", "GIR 50% 图", "DZiG/GraphBolt 50% 图"],
       "BC": ["GIR", "GIR 对称图", "Terrace", "Aspen", "Spruce"]}
for a in ("BFS", "SSSP", "BC"):
    w(f"### 2.{'BFS SSSP BC'.split().index(a) + 1} {a}\n")
    w("| 系统 | " + " | ".join(ABBR) + " |")
    w("|---|" + "---|" * len(ABBR))
    for s in SYS[a]:
        w(f"| {s} | " + " | ".join(f(alg(ds, s, a)) for ds, _ in DATASETS) + " |")
    w("")
# speedups vs paper
PAIRS = [("BFS", "Terrace", "GIR", "Terrace"), ("BFS", "Aspen", "GIR 对称图", "Aspen"), ("BFS", "Spruce", "GIR", "Spruce"),
         ("BFS", "GraphBolt", "GIR 50% 图", "DZiG/GraphBolt 50% 图"), ("BFS", "DZiG", "GIR 50% 图", "DZiG/GraphBolt 50% 图"),
         ("SSSP", "Spruce", "GIR", "Spruce"), ("SSSP", "GraphBolt", "GIR 50% 图", "DZiG/GraphBolt 50% 图"),
         ("SSSP", "DZiG", "GIR 50% 图", "DZiG/GraphBolt 50% 图"),
         ("BC", "Terrace", "GIR 对称图", "Terrace"), ("BC", "Aspen", "GIR 对称图", "Aspen")]
w("### 2.4 GIR 加速比：论文表 4 与实测\n")
w("对方时间 ÷ GIR 时间（大于 1 表示 GIR 更快）。Aspen、Terrace 的 BC 与对称图上的 GIR 比；DZiG/GraphBolt 与同一 50% 图上的 GIR 从头计算比。\n")
w("| 算法 | 对比 | 来源 | " + " | ".join(ABBR) + " | 平均 |")
w("|---|---|---|" + "---|" * (len(ABBR) + 1))
for a, name, gsys, osys in PAIRS:
    vals, cells = [], []
    for ds, _ in DATASETS:
        g, o = alg(ds, gsys, a), alg(ds, osys, a)
        cells.append(ratio(o, g))
        vals.append(o / g if g and o and g > 0 and o > 0 else None)
    gm = geo(vals)
    pv = P_ALG[(a, name)]
    w(f"| {a} | {name} | 论文 | " + " | ".join(f"{x}×" for x in pv[:7]) + f" | {pv[7]}× |")
    w(f"| {a} | {name} | 实测 | " + " | ".join(cells) + f" | {ratio(gm, 1.0) if gm else '—'} |")
w("\nGraphBolt 与 DZiG 的 BFS/SSSP 是同一个引擎，所以两行实测相同；论文给了不同的数字。\n")

# ---------------------------------------------------------------- scalability
w("## 3. BFS 线程扩展性（论文第 7.4 节，图 9）\n")
w("各系统相对自身单线程的加速比（单线程时间 ÷ N 线程时间）。TW、FS 只测 1、2、8、32、64 线程。\n")
SC = [("GIR", lambda ds: (ds, "GIR", "bfs_ms", "direction-optimizing")),
      ("Terrace", lambda ds: (ds, "Terrace", "bfs_ms", "push-only")),
      ("Aspen", lambda ds: (ds + "_sym", "Aspen", "bfs_ms", "stock")),
      ("Spruce", lambda ds: (ds, "Spruce", "bfs_ms", "top-down"))]
w("| 系统 | 数据集 | " + " | ".join(f"{t} 线程" for t in (1, 2, 4, 8, 16, 32, 64)) + " | 64 线程加速比 | 论文 |")
w("|---|---|" + "---|" * 9)
avg = {}
for sname, key in SC:
    for ds, ab in DATASETS:
        k = key(ds)
        ts = {t: get(*k, threads=t) for t in (1, 2, 4, 8, 16, 32, 64)}
        if all(v is None for v in ts.values()):
            continue
        t1 = ts.get(1)
        ok = lambda v: v is not None and v > 0
        cells = [f(ts[t]) + (f" ({t1 / ts[t]:.1f}×)" if ok(t1) and ok(ts[t]) and t != 1 else "") if ts[t] is not None else "—" for t in (1, 2, 4, 8, 16, 32, 64)]
        sp = t1 / ts[64] if ok(t1) and ok(ts.get(64)) else None
        avg.setdefault(sname, []).append(sp)
        pc = P_SCALE.get(sname, {}).get(ab)
        w(f"| {sname} | {ab} | " + " | ".join(cells) + f" | {f'{sp:.1f}×' if sp else '—'} | {f'{pc}×' if pc else '—'} |")
w("\n格中为毫秒，括号内为相对单线程的加速比。\n")
w("| 系统 | 64 线程平均加速比（实测，几何平均） | 论文 |")
w("|---|---|---|")
for sname, _ in SC:
    gm = geo(avg.get(sname, []))
    pc = P_SCALE.get(sname, {}).get("Avg")
    w(f"| {sname} | {f'{gm:.1f}×' if gm else '—'} | {f'{pc}×' if pc else '—（论文只给 32 线程：Terrace 15.5×，Aspen 15.7×）'} |")
w("")

# ---------------------------------------------------------------- streaming
w("## 4. 流式场景（DZiG 的设定：每批先改图，再求该源点的 BFS/SSSP）\n")
w("每批合计时间（ms，越低越好），10 批 × 10,000 次操作，源点相同。GIR 从头计算与增量修正两种；“延迟维护”每批只写边存储。\n")
w("| 数据集 | 算法 | 流 | DZiG | GIR 批量 + 从头 | GIR 批量 + 增量 | GIR 延迟 + 增量 | DZiG ÷ GIR 延迟 + 增量 |")
w("|---|---|---|---|---|---|---|---|")
for ds, ab in DATASETS:
    for a in ("bfs", "sssp"):
        for s in ("insert", "delete", "mixed"):
            m = f"{a}_{s}_endtoend_ms"
            dz = get(ds, "DZiG", m)
            g0 = get(ds, "GIR", m)
            gi = get(ds, "GIR-inc", m)
            gl = get(ds, "GIR-lazy-inc", m)
            if dz is None and g0 is None:
                continue
            w(f"| {ab} | {a.upper()} | {s} | {f(dz)} | {f(g0)} | {f(gi)} | {f(gl)} | {ratio(dz, gl)} |")
w("")
mm = [(ds, get(ds, s, f"inc_mismatch_{st}")) for ds, _ in DATASETS for s in ("GIR-inc", "GIR-lazy-inc") for st in ("insert", "delete", "mixed")]
bad = [x for x in mm if x[1]]
w(f"增量结果在每条流结束时与从头计算逐点比对：{'全部一致' if not bad else '有不一致：' + str(bad)}。\n")

# ---------------------------------------------------------------- correctness
w("## 5. 结果一致性\n")
w("第一个源点的 BFS 可达顶点数与 SSSP 距离和（同一图上的系统应相同）。\n")
w("| 数据集 | GIR 可达 | Terrace 可达 | Spruce 可达 | GIR 距离和 | Spruce 距离和 | 对称图：GIR / Terrace / Aspen 可达 |")
w("|---|---|---|---|---|---|---|")
for ds, ab in DATASETS:
    def g(dd, s, m):
        v = get(dd, s, m)
        return f"{v:,.0f}" if v is not None and not math.isnan(v) else "—"
    w(f"| {ab} | {g(ds, 'GIR', 'bfs_reached_src0')} | {g(ds, 'Terrace', 'bfs_reached_src0')} | {g(ds, 'Spruce', 'bfs_reached_src0')} | "
      f"{g(ds, 'GIR', 'sssp_distsum_src0')} | {g(ds, 'Spruce', 'sssp_distsum_src0')} | "
      f"{g(ds + '_sym', 'GIR', 'bfs_reached_src0')} / {g(ds + '_sym', 'Terrace', 'bfs_reached_src0')} / {g(ds + '_sym', 'Aspen', 'bfs_reached_src0')} |")
w("")

lp = os.path.join(PD, "load.csv")
if os.path.exists(lp):
    w("## 6. 运行时的机器负载\n")
    w("每一步开始前等待机器空闲（1 分钟负载 < 4，其他用户 CPU < 1.5 核，连续 3 次检查），运行中每 10 秒采样其他用户的 CPU 占用；"
      "超过 3 核的步骤等空闲后重跑一次。下表是各数据集的汇总。\n")
    st = {}
    for r in csv.DictReader(open(lp)):
        s = st.setdefault(r["dataset"], [0, 0, 0, 0.0, 0])
        s[0] += 1; s[1] += int(r["disturbed"]); s[2] = max(s[2], int(r["max_other_cpu_pct"])); s[3] += float(r["seconds"]); s[4] += int(r["attempt"]) == 2
    w("| 数据集 | 步骤数（含重跑） | 受扰 | 其中重跑 | 其他用户 CPU 最大值（%） | 总时长（h） |")
    w("|---|---|---|---|---|---|")
    for ds, s in st.items():
        w(f"| {ds} | {s[0]} | {s[1]} | {s[4]} | {s[2]} | {s[3] / 3600:.1f} |")
    w("")
open(OUTF, "w").write("\n".join(out) + "\n")
print("wrote", OUTF)
