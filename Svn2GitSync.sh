#!/bin/bash

INITIALIZE=$1
SVNREPOBASENAME=$2
IGNORES=$3
SVNWORKSPACEPATH=$HOME/workspace/$SVNREPOBASENAME
TMPGITPATH=$HOME/tmp/git

echo Usage:
echo $0 "\"init/rebase\" \"svn repo basename\" \"ignores (comma separated!)\""

if [ $# -eq 0 ]; then
    echo "Wrong number of paramers! => EXITING..."
    exit
fi

if [ "$INITIALIZE" == "init" ];then
    Arguments="-s svn+ssh://clio.psy.vu.nl/home/r.pool/svn/$SVNREPOBASENAME $SVNREPOBASENAME --username=r.pool"
    if [ "$IIGNORES"!="" ];then
	    for i in `echo $IGNORES | sed -e s/","/" "/g`
	    do
	        Arguments=$Arguments" --ignore-paths=$i"
	    done
	fi
	cd $TMPGITPATH
	if [ "$SVNREPOBASENAME"!="" ]; then
	    rm -rf $SVNREPOBASENAME
	fi
    git svn clone -s $Arguments
    cd $SVNREPOBASENAME
    git remote add origin git@github.com:rpool/$SVNREPOBASENAME.git
    git push origin --all
    git push origin --tags
    for i in `echo $IGNORES | sed -e s/","/" "/g`
	do
	    echo $i >> .git/info/exclude
	    echo .svn >> .git/info/exclude
	done
    rm -rf $SVNWORKSPACEPATH/.git
    cp -rp .git $SVNWORKSPACEPATH
fi

if [ "$INITIALIZE" == "rebase" ]; then
    cd $SVNWORKSPACEPATH
    if [ "$IIGNORES"!="" ];then
        for i in `echo $IGNORES | sed -e s/","/" "/g`
        do
            Arguments=$Arguments" --ignore-paths=$i"
        done
    fi
    git reset --hard HEAD
    git svn rebase $Arguments
    git push origin --all
    git push origin --all
fi

echo "All done ..."
echo "Bye bye ;-)"

exit 0