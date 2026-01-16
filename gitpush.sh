#!/bin/bash

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# 显示菜单
show_menu() {
    echo -e "${BLUE}==============================${NC}"
    echo -e "         Git 助手脚本         "
    echo -e "${BLUE}==============================${NC}"
    echo "当前本地分支: ${GREEN}$BRANCH${NC}"
    echo "1. 查看状态 (git status)"
    echo "2. 拉取更新 (git pull)"
    echo "3. 提交并推送 (add + commit + push)"
    echo "4. 仅推送到指定远程分支"
    echo "0. 退出"
    echo -e "${BLUE}==============================${NC}"
}

# 获取当前分支
BRANCH=$(git symbolic-ref --short HEAD)

while true; do
    show_menu
    read -p "请输入选项 [0-4]: " choice
    
    case $choice in
        1)
            git status
            ;;
        2)
            echo -e "${YELLOW}正在从 origin/${BRANCH} 拉取代码...${NC}"
            git pull origin "$BRANCH"
            ;;
        3)
            git status -s
            echo ""
            read -p "请输入 Commit 信息 (直接回车默认 'Update'): " msg
            msg=${msg:-Update}

            read -p "请输入要推送的远程分支名 (直接回车默认 $BRANCH): " target_branch
            target_branch=${target_branch:-$BRANCH}

            echo -e "${GREEN}执行: git add .${NC}"
            git add .
            
            echo -e "${GREEN}执行: git commit -m \"$msg\"${NC}"
            git commit -m "$msg"
            
            echo -e "${GREEN}执行: git push origin $BRANCH:$target_branch${NC}"
            git push origin "$BRANCH:$target_branch"
            ;;
        4)
            read -p "请输入要推送的远程分支名 (直接回车默认 $BRANCH): " target_branch
            target_branch=${target_branch:-$BRANCH}
            echo -e "${GREEN}执行: git push origin $BRANCH:$target_branch${NC}"
            git push origin "$BRANCH:$target_branch"
            ;;
        0)
            echo "Bye!"
            exit 0
            ;;
        *)
            echo -e "${RED}无效选项，请重试${NC}"
            ;;
    esac
    echo ""
    read -p "按回车键继续..."
done
