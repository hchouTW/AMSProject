"/**************************************************************
" *  Author        : Hsin-Yi Chou
" *  Email         : hchou@cern.ch
" *  Version       : 1.0
" *  Last modified : 28 Mar 2017 
" *  Filename      : .vimrc
" *  Description   : VIM setting
" *************************************************************/


" Reference Author : Doug Black
" https://dougblack.io/words/a-good-vimrc.html 


" Color Scheme (Solarized)
" http://ethanschoonover.com/solarized
" https://github.com/altercation/Vim-colors-solarized
" move solarized.vim to ~/.vim/colors/
if has ('gui_running')
	set background=light
else
	set background=light
endif
colorscheme solarized

set t_Co=256
syntax on " enable syntas processing


" Spaces & Tabs
set tabstop=4     " number of visual spaces per TAB
set softtabstop=4 " number of spaces in tab when editing
set shiftwidth=4  " number of spaces when indenting
"set expandtab     " tabs are spaces 


" UI Config
set number     " show line numbers
set showcmd    " show command in bottom bar
set cursorline " highlight current line
set wildmenu   " visual autocomplete
set showmatch  " highlight matching {[()]}

filetype indent on " load filetype-specific files


" Searching
set incsearch  " search as characters are entered
set hlsearch   " highlight matches


" Folding
set foldenable        " enable folding
set foldlevelstart=10 " open most folds by default
set foldnestmax=10    " 10 nested fold max
set foldmethod=indent " fold based on indent level

" space open/close folds
nnoremap <space> za


" Movement
" move vertically by visual line
nnoremap j gj
nnoremap k gk

" move to beginning/end of line
nnoremap B ^
nnoremap E $

" $/^ doesn't do anything
nnoremap ^ <nop>
nnoremap $ <nop>


" Mouse
"set mouse=nv


" Backups
set backup


" in normal mode F2 will save the file
nmap <F2> :w<CR>
" in insert mode F2 will exit insert, save, enters insert again
imap <F2> <ESC>:w<CR>i
" switch between header/source with F4
map <F4> :e %:p:s,.h$,.X123X,:s,.C$,.h,:s,.X123X$,.C,<CR>

" Enhanced keyboard mappings
nmap <f4> :e %:r.h<cr>
