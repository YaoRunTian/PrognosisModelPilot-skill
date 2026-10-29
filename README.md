# PrognosisModelPilot-skill

面向肿瘤预后研究的 Codex skill，整合单细胞与转录组数据，支持项目方案、R 分析项目、结果图表和实验报告。

技能入口位于 [prognosis-model-pilot/SKILL.md](prognosis-model-pilot/SKILL.md)。完整项目默认从方案和技术路线开始，再生成 R 代码、执行分析并依据实际结果撰写报告。

研究主线包括关键细胞基因集、转录组基因集与可选通路基因集的整合，单因素 Cox、LASSO-Cox 和多因素 Cox 建模，以及模型验证、临床列线图和相应下游分析。数据缺失或分析未完成时，交付物需标明实际状态。

## 目录

- `prognosis-model-pilot/SKILL.md`：技能说明与任务入口。
- `prognosis-model-pilot/references/`：流程、部门标准及参考项目整理。
- `prognosis-model-pilot/assets/`：方案、报告和交付物清单模板。
- `prognosis-model-pilot/scripts/`：TCGA 数据获取和交付文件检查脚本。

本仓库附带五个项目的根目录R/Rmd参考脚本；未包含完整服务器代码库或患者数据。公开版隐藏本地与服务器绝对路径，不包含部门标准原件；标准提炼仍在 references/ 中，原件需由使用者自行提供。参考脚本中的路径占位符及项目专用参数需要适配；参考代码经过语法检查，尚未作为完整流程在真实数据上复现。
