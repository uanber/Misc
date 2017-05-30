# ~/.cshrc for unified home directory (USED FOR GFDL HOME DIRECTORY)

source /home/gfdl/init/csh.cshrc

# if remote command or batch job, avoid aliases below
if ( "`tty`" == "not a tty" ) exit

#-------------DO NOT CHANGE ABOVE THIS LINE-------------

# script run from command line
if ( ! $?prompt ) exit


# login shell, interactive sub-shell

# tcsh filename completion via TAB key
# unset addsuffix   # turn off trailing slash or space
# unset autolist    # turn off listing of alternatives (use ^D to list)


# settings for generic Linux (workstations or HPCS)
alias ls  'ls -F'
alias du  'du -h'
alias df  'df -h'

if ( `gfdl_platform` == desktop ) then
  # settings for workstations only
endif

if ( `gfdl_platform` == hpcs-csc ) then
  # settings for HPCS only
endif

stty erase '^?'
