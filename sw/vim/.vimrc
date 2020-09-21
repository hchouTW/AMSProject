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


" Backups
" save backup file in backupdir=./.bck
let g:locbakdir=resolve(expand('%:p:h') . "/.bak")
autocmd BufEnter    * let g:locbakdir=resolve(expand('%:p:h') . "/.bak")
autocmd BufWritePre * execute '!mkdir -p ' . g:locbakdir 
autocmd BufWritePre * :let &backupdir=g:locbakdir
set backup


" Color Scheme (Solarized)
" http://ethanschoonover.com/solarized
" https://github.com/altercation/Vim-colors-solarized
" move solarized.vim to ~/.vim/colors/
if has ('gui_running')
    set background=light
else
    set background=light
endif
let g:solarized_termcolors=256
" colorscheme solarized
syntax on " enable syntas processing

""" vim 8.0
" packadd! dracula
" syntax enable
" colorscheme dracula

"" 120 Column Layout
"set textwidth=120
"if exists('+colorcolumn')
"    set colorcolumn=120
"else
"    au BufWinEnter * let w:m2=matchadd('ErrorMsg', '\%>120v.\+', -1)
"endif


" Spaces & Tabs
set tabstop=4     " number of visual spaces per TAB
set softtabstop=4 " number of spaces in tab when editing
set shiftwidth=4  " number of spaces when indenting
set expandtab     " tabs are spaces 
set autoindent    " copy indent from current line when starting a new line 
set smartindent   " do smart autoindenting when starting a new line


" UI Config
set number      " show line numbers
set showcmd     " show command in bottom bar
set cursorline  " highlight current line
set wildmenu    " visual autocomplete
set showmatch   " highlight matching {[()]}
set scrolloff=5 " minimal number of screen lines to keep above and below the cursor

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

" in normal mode F8 will use internal formatting
nmap <F8> ggVG=

" in normal mode F2 will save the file
nmap <F2> :w<CR>

" in insert mode F2 will exit insert, save, enters insert again
imap <F2> <ESC>:w<CR>i

" switch between header/source with F4
map <F4> :e %:p:s,.h$,.X123X,:s,.C$,.h,:s,.X123X$,.C,<CR>

" Enhanced keyboard mappings
nmap <f4> :e %:r.h<cr>



let g:tex_flavor='latex'
