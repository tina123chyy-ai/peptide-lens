# Peptide Lens · 突变肽序列分析工具

**C++17 · Data structures · Bioinformatics / C++17 · 数据结构 · 生物信息**

A small, reproducible command-line project for comparing paired protein sequences and extracting mutation-containing peptide windows. / 一个可复现的小型命令行项目：比较配对的蛋白质序列，提取包含变化位置的短肽片段。

> **Scope / 定位**: Educational sequence preprocessing, not a validated neoantigen predictor. All bundled sequences and substitutions are synthetic. No patient data, wet-lab results, HLA predictions, vaccine efficacy predictions or CAR-T efficacy predictions are included. / 本项目用于教学和序列预处理，不是经过验证的新抗原预测工具。附带序列与突变全部为模拟数据，不含患者数据、实验结果、HLA 预测或疫苗及 CAR-T 疗效预测。

## 1. What it does / 它做什么

Give the program two FASTA files: a normal sequence file and a corresponding mutant sequence file. It matches records by ID, identifies amino-acid substitutions, extracts windows covering those substitutions, groups repeated mutant peptides and writes CSV tables. / 输入正常与突变两个 FASTA 文件；程序按编号配对，找出氨基酸替换，截取覆盖替换位置的肽段，合并重复肽段并导出 CSV 表格。

This connects to the sequence-preprocessing stage of candidate neoantigen research. It does not implement the combination therapy described in the associated review topic. / 本项目与候选新抗原研究中的序列预处理有关，不实现相关综述中讨论的联合治疗。

## 2. Mac quick start / Mac 快速开始

**You do not need GitHub to run this project. / 不需要 GitHub 账号也能运行。**

1. Unzip the package. / 解压下载包。
2. Open Terminal using Spotlight (Command + Space, type “Terminal”). / 用聚焦搜索打开“终端”。
3. If C++ developer tools are not installed, run the following, complete the installation dialog and then continue. / 如果没有安装 C++ 开发工具，运行下面命令，完成弹窗中的安装后继续：

```bash
xcode-select --install
```

If it says the tools are already installed, continue. Full Xcode and Homebrew are not required. / 若提示已安装，继续即可；无需完整 Xcode，也无需 Homebrew。

4. Type `cd ` (including the trailing space), drag the extracted `peptide-lens` folder into Terminal, then press Enter. / 在终端输入 `cd `（末尾有空格），将解压后的 `peptide-lens` 文件夹拖进去，再按回车。
5. Run / 输入：

```bash
bash run_demo.command
```

The script compiles locally, runs the demo and opens the results folder on Mac. The same source is intended for Apple Silicon and Intel Macs through local compilation. / 脚本在本机编译并运行示例，完成后在 Mac 上打开结果文件夹。同一份源码通过本机编译适配 Apple 芯片和 Intel Mac。

**Verification boundary / 验证范围**: The code and tests were executed in Linux with GCC, not on a physical Mac. The Mac launcher uses Bash and standard C++17 compilation; actual macOS execution remains to be confirmed on the recipient's computer. No Linux executable is shipped as a Mac binary. / 本次已在 Linux/GCC 环境执行代码和测试，未在实体 Mac 上实测。Mac 启动脚本采用 Bash 和标准 C++17 编译；仍需在接收者的 Mac 上确认运行。不将 Linux 可执行文件冒充 Mac 软件。

If double-clicking the `.command` file is blocked, use the Terminal command above; do not disable system security. / 如双击脚本无法打开，使用上面的终端命令即可，不需要关闭系统安全设置。

## 3. Example results / 示例结果

**Start with `results/00_结果阅读说明.txt` for plain-language file descriptions and column meanings. / 建议先打开 `results/00_结果阅读说明.txt`，查看每个文件回答什么问题及各列含义。**

Default peptide lengths are 8, 9, 10 and 11. They are configurable extraction lengths, not biological eligibility criteria. / 默认提取长度为 8、9、10、11；这只是可配置的片段长度，不是生物学筛选或有效性判断标准。

The bundled dataset contains five pairs of 40-amino-acid sequences: / 示例含五对长度为 40 的序列：

| ID | Setup / 设置 |
|---|---|
| demo_single | One substitution at 15 / 第15位单点替换 |
| demo_boundary | Substitutions at 1 and 40 / 首尾替换 |
| demo_multi | Substitutions at 15 and 18 / 第15及18位替换 |
| demo_control | No substitutions / 无变化对照 |
| demo_duplicate | Same sequence pair as demo_single / 与单点样本相同，用于验证去重和来源保留 |

**Verified demo: 5 samples, 6 substitutions, 134 candidate occurrences, 84 unique mutant peptides; the sample query returns 38 rows. / 实测示例：5个样本、6处替换、134条片段记录、84条不同突变肽；按单点样本查询得到38条记录。**

See `results/demo_summary.txt` for actual counts and `data/expected_mutations.csv` for the known answer. / 实际数量见 `results/demo_summary.txt`；事先设定的答案见 `data/expected_mutations.csv`。

| Output / 文件 | Meaning / 含义 |
|---|---|
| demo_mutations.csv | One row per substitution / 每行一个替换位置 |
| demo_candidates.csv | All mutation-containing windows, with normal counterpart and source / 所有含替换的片段，保留对应正常肽和来源 |
| demo_unique.csv | Identical mutant peptides grouped, with occurrence and sample counts / 相同突变肽合并，列出出现次数及样本数 |
| demo_queries.csv | Matching rows from an optional sample or peptide query / 按样本或肽序列查询所得记录 |
| demo_summary.txt | Counts and extraction settings / 数量及提取设置 |
| test_report.txt | Actual automated-test log / 实际自动测试记录 |

All coordinates are **1-based and inclusive**. Mutation positions separated by `;` refer to full-sequence coordinates, not relative peptide positions. / 坐标从 **1** 开始，包含首尾；用分号分隔的突变位置指整条蛋白序列中的位置，不是片段内的相对位置。

An occurrence is a sample + start + length combination. A window covering multiple mutations appears only once, with all changes listed. Deduplication never deletes occurrence records. / 一条片段记录对应“样本＋起点＋长度”；覆盖多个突变的同一窗口只记录一次，并列出全部变化。去重不删除来源明细。

## 4. Custom runs / 自定义运行

Build once / 先编译：

```bash
bash build.sh
```

Run with alternative lengths / 改变提取长度：

```bash
./build/peptide-lens --normal data/normal.fasta --mutant data/mutant.fasta --lengths 8,9 --out results/custom
```

Query one sample / 按样本查询：

```bash
./build/peptide-lens --normal data/normal.fasta --mutant data/mutant.fasta --out results/sample --sample demo_single
```

Query an exact, uppercase mutant peptide / 精确查询大写突变肽：

```bash
./build/peptide-lens --normal data/normal.fasta --mutant data/mutant.fasta --out results/search --query IKLMNPQA
```

`--sample` and `--query` may be combined (intersection). They filter only `_queries.csv`, not the full analysis. Without filters, `_queries.csv` contains a header only. No matches are a valid zero-result query. / 两个查询条件可同时使用（取交集），只影响查询表，不改变全量结果；不指定查询条件时查询表仅有表头。无匹配是合法的零结果。

Output folders must already exist; an existing identical output prefix is overwritten. Use a new prefix to retain old runs. / 输出文件夹须已存在；相同前缀会覆盖已有结果，保留旧结果请换前缀。

## 5. Input format / 输入格式

FASTA headers begin with `>`; the first whitespace-delimited token is the ID. The two files must contain exactly the same unique IDs; order may differ. Each pair must have equal length and be positionally aligned. / FASTA 标题以 `>` 开头，第一个空白分隔的字符串作为 ID。两个文件必须包含完全相同且不重复的 ID，顺序可以不同。每一对须等长且位置已经对应。

Only the 20 standard amino-acid letters `ACDEFGHIKLMNPQRSTVWY` are accepted. Lowercase is converted to uppercase; multiline sequences are accepted. Unknown residues, gaps, stop symbols, empty sequences and unequal lengths are rejected. / 接受20种标准氨基酸字母，小写自动转为大写，可分多行；不接受未知残基、缺口、终止符、空序列及不等长序列。

This is a positional comparison, **not a sequence-alignment algorithm**. Insertions, deletions, frameshifts, DNA/RNA input and ambiguous residues are outside version 1. / 本版本是逐位比较，**不执行序列比对校准算法**；不支持插入、删除、移码、DNA/RNA 输入或模糊残基。

## 6. C++ and data structures / C++ 与数据结构

| Component / 组件 | Purpose / 用途 |
|---|---|
| struct | Sequence, Mutation, Candidate records / 序列、替换、候选片段记录 |
| vector | Ordered lists of records and mutation positions / 记录列表和变化位置 |
| unordered_map | Match sample IDs and index identical peptides / 样本匹配、肽序列哈希索引 |
| set | Validate unique IDs, lengths and sample counts / 唯一编号、长度及样本计数 |
| Sliding windows + lower_bound | Find windows overlapping substitutions / 查找覆盖替换位置的窗口 |
| sort | Stable, deterministic output order / 确定性输出顺序 |

Per pair of length N with M substitutions, scanning costs O(N); window selection across K lengths costs O(K N log(M+1)), plus the cost of copying emitted peptides and their mutation positions. Hash lookup is expected O(1) in entry count, with string hashing proportional to peptide length; worst-case collisions can degrade it. / 每对长度 N、含 M 个变化的序列，逐位扫描为 O(N)；K 种长度的窗口筛选约为 O(K N log(M+1))，另加输出肽及突变位置复制开销。哈希查询按条目数平均为 O(1)，但字符串哈希仍与肽长有关，碰撞最坏情况会退化。

## 7. Testing / 测试

The demo needs only a C++ compiler. Optional automated tests additionally require Python 3, with no third-party packages. / 示例仅需 C++ 编译器；可选自动测试另外需要 Python 3，不需要第三方包。

```bash
bash build.sh
python3 tests/test_cli.py
```

Tests cover known mutation answers, sequence ends, no-change controls, overlapping mutations, duplicate peptides, sample matching, CSV escaping, invalid inputs and 80 deterministic randomized comparisons against an independent window-enumeration oracle. / 测试覆盖已知答案、首尾、无变化、多突变重叠、重复肽、样本配对、CSV 转义、错误输入及80组固定随机种子的独立窗口枚举核对。

## 8. Limitations and provenance / 局限与来源

- All sample sequences are fabricated by repeating the standard amino-acid alphabet; they are not database proteins or patient observations. See `data/PROVENANCE.md`. / 所有示例由标准氨基酸字母表重复构造，不是数据库蛋白或患者观测，详见数据说明。
- A mutation-containing peptide is only a sequence candidate. There is no HLA binding, expression, processing, self-proteome filtering, immune recognition, clinical safety or efficacy assessment. / 含突变的肽仅是序列层面的候选；不评估 HLA 结合、表达、加工、自身蛋白组匹配、免疫识别或临床安全有效性。
- No claim of original experimental findings or independent student authorship is made. The first implementation and documentation were generated with AI assistance and executed against synthetic tests. Describe actual contributions honestly when presenting this project. / 不声称有原创实验发现或学生独立完成经历。首版代码与文档由 AI 辅助生成并用模拟数据测试；展示时应如实说明实际贡献。

## 9. GitHub handoff / GitHub 交付

Upload the extracted project contents (not `build/`) to the repository chosen later. `README.md` will be the bilingual landing page. No account, credentials or remote repository is bundled. / 后续将解压后的项目内容（不含 build 目录）上传到指定仓库；README.md 自动作为双语首页。本包不包含账号、密码或远程仓库。
