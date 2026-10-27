# Data provenance / 数据来源

All sequences are **entirely synthetic**, created for software verification. No external protein sequence database was used. / 全部序列均为软件验证而人工构造，未使用外部蛋白序列数据库。

The normal template is `ACDEFGHIKLMNPQRSTVWY` repeated twice (40 residues). / 正常模板是20种标准氨基酸字母表重复两次，共40位。

Substitutions use 1-based positions: / 替换位置从1开始计数：

- demo_single / demo_duplicate: R15A.
- demo_boundary: A1Y, Y40A.
- demo_multi: R15A, V18C.
- demo_control: unchanged / 无变化。

`expected_mutations.csv` was prepared from these specified edits before testing. The deliberately repeated alphabet is not biologically representative, so counts, repeat frequencies and performance on this tiny dataset must not be interpreted as tumor biology. / 答案文件按事先指定的变化建立。重复字母模板不具备生物学代表性，不能将其数量或重复频率解释成真实肿瘤特征。

FASTA headers also label each record synthetic. / FASTA 标题同样注明模拟性质。
