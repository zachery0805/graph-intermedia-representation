#!/bin/bash
# 自动上传到 GitHub 的脚本

# 获取当前时间作为提交信息
TIME=$(date "+%Y-%m-%d %H:%M:%S")

# 添加所有修改
git add .

# 提交（自动加时间）
git commit -m "Auto update at $TIME"

# 拉取远程更新，避免冲突
git pull origin main --allow-unrelated-histories

# 推送到远程仓库
git push origin main
