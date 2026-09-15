# Experiment records

每个数值实验必须先得到用户对精确参数的确认，随后在服务器的 `runs/<experiment-id>/` 下建立
独立记录。`runs/` 不进入 Git，但记录必须在 AutoDL-TMP 工作期内保持完整。

每个实验目录至少包含：

- `manifest.yaml`：目的、结论边界、代码提交、环境、输入路径、参数、并行度和时间；
- `input_sha256.txt`：实际读取的输入文件及其 SHA-256 校验值；
- `stdout.log`：未经删节的标准输出与错误输出；
- `raw/`：求解器直接产物；
- `figures/`：由该实验数据生成的检查图；
- `summary.md`：是否成功、实际运行范围、异常和对论文图的影响。

记录原则：失败或中止的实验也保留；不得覆盖旧实验目录；未经新确认不得扩大扫描范围；正式论文图
只从 manifest 指向的原始数据生成。

实验完成后，将便于审阅的精简副本发布到 `experiments/results/<experiment-id>/` 并同步当前
工作分支到 GitHub。精简副本保留 README、manifest、核心数值 CSV、对照图与必要的校验信息；
体积较大的完整运行日志只留在 AutoDL-TMP 的 `runs/<experiment-id>/`。
