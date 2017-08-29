# AMSProject
AMS02 Project

## Fast Setting
In ~/.bashrc file
```bash
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
export AMSProj=~/AMSProject
export AMSCore=~/AMSCore

source ${AMSProj}/sys/shell/bashrc.sh
```
VIM
```bash
ln -s ${AMSProj}/sw/vim ~/.vim
ln -s ${AMSProj}/sw/vim/.vimrc ~/.vimrc
```
ROOT
```bash
ln -s ${AMSProj}/sw/ROOT/.rootlogon.C .rootlogon.C
```
Git
```bash
ln -s ${AMSProj}/sw/git/.gitconfig .gitconfig
```

## ENV Setting
read env/README.md

## VIM Setting
read vim/README.md

## Coding Style
read guides/READE.md
