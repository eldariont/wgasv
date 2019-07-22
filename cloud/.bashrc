#Aliases
alias ..='cd ..'
alias ...='cd ../..'
alias la='ls -lisah'
alias gl='git log --graph --all --pretty=format:"%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%r) %C(bold blue)<%an>%Creset" --abbrev-commit --date=relative'
alias hg='history | grep'

umask 002

#Set PS1 string
if [ "$TERM" = "linux" ]
then
    #we're on the system console or maybe telnetting in
    export PS1="\[\e[32;1m\]\u@\H > \[\e[0m\]"
else
    #we're not on the console, assume an xterm
    export PS1="\[\e]2;\u@\H \w\a\e[32;1m\]\w>\[\e[0m\] " 
fi

export HISTCONTROL=ignoredups:erasedups  # no duplicate entries
export HISTSIZE=100000                   # big big history
export HISTFILESIZE=100000               # big big history
shopt -s histappend                      # append to history, don't overwrite it

#PATH
PATH="/ebs1/vcflib/bin:$PATH"
PATH="/ebs1/seqwish/bin:$PATH"
PATH="/ebs1/minimap2:$PATH"
PATH="/ebs1/minimap2/misc:$PATH"
PATH="/ebs1/bcftools:$PATH"
PATH="/ebs1/hal2vg:$PATH"
PATH="/ebs1/repeatMaskerPipeline:$PATH"
PATH="/ebs1/AsmVar/src/AsmvarDetect:$PATH"
PATH="/ebs1/vg/bin:/ebs1/vg/scripts:$PATH"
export PATH
export LIBRARY_PATH=`pwd`/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=`pwd`/lib:$DYLD_LIBRARY_PATH
export LD_INCLUDE_PATH=`pwd`/include:$LD_INCLUDE_PATH
export C_INCLUDE_PATH=`pwd`/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=`pwd`/include:$CPLUS_INCLUDE_PATH
export INCLUDE_PATH=`pwd`/include:$INCLUDE_PATH

export CC=$(which gcc)
export CXX=$(which g++)
