#!/bin/bash
# =========================================================
# 自动上传到 GitHub 的脚本（带认证和错误处理）
# 作者: Chaoqun Li
# 日期: $(date "+%Y-%m-%d %H:%M:%S")
# =========================================================

set -e  # 如果有错误立即退出
TIME=$(date "+%Y-%m-%d %H:%M:%S")

echo "🔍 检查 Git 远程仓库连接方式..."
REMOTE_URL=$(git remote get-url origin)

# --- 1. 检查是否使用 SSH ---
if [[ $REMOTE_URL == git@github.com:* ]]; then
    echo "✅ 检测到 SSH 连接方式，将直接推送。"
else
    echo "⚠️ 远程使用 HTTPS，将尝试使用凭证助手。"

    # 确保凭证助手可用
    git config --global credential.helper store

    # 检查是否已有凭证缓存
    if ! git credential fill <<EOF | grep -q 'password='
protocol=https
host=github.com
EOF
    then
        echo "🔑 未检测到已缓存的 GitHub 凭证。"
        echo "请输入你的 GitHub 用户名（例如：Kepavel）："
        read USERNAME
        echo "请输入你的 GitHub Personal Access Token（PAT）："
        read -s TOKEN
        git remote set-url origin "https://${USERNAME}:${TOKEN}@github.com/Kepavel/gir.git"
        echo "✅ 已设置临时凭证。"
    fi
fi

# --- 2. 添加、提交、拉取、推送 ---
echo "📦 添加修改..."
git add .

echo "📝 提交更新..."
git commit -m "Auto update at $TIME" || echo "（无改动可提交）"

echo "⬇️ 拉取远程更新以避免冲突..."
git pull origin main --allow-unrelated-histories --rebase

echo "⬆️ 推送到 GitHub..."
git push origin main

echo "🎉 推送完成于 $TIME"
