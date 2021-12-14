#!/bin/bash

NAME='RTCpredictor'
DIR=~/.$NAME
if [ -d $DIR ]
then
	rm -rf $DIR
	mkdir $DIR
else
	mkdir $DIR
fi
wget https://cpan.metacpan.org/authors/id/E/ET/ETHER/Try-Tiny-0.30.tar.gz
tar -xzf Try-Tiny-0.30.tar.gz
cd Try-Tiny-0.30
perl Makefile.PL PREFIX=$DIR
make
make install
cd ../

LDIR=`find $DIR -name Tiny.pm`
LDIR=${LDIR/\/Try\/Tiny.pm/}
export PERL5LIB=$LDIR:$DIR/lib64/perl5:$DIR/share/perl5

wget https://cpan.metacpan.org/authors/id/R/RJ/RJBS/Test-Fatal-0.016.tar.gz
tar -xzf Test-Fatal-0.016.tar.gz
cd Test-Fatal-0.016
perl Makefile.PL PREFIX=$DIR
make
make install
cd ../

wget https://cpan.metacpan.org/authors/id/H/HA/HAARG/Test-Needs-0.002009.tar.gz
tar -xzf Test-Needs-0.002009.tar.gz
cd Test-Needs-0.002009
perl Makefile.PL PREFIX=$DIR
make
make install
cd ../

wget https://cpan.metacpan.org/authors/id/E/ET/ETHER/Class-Method-Modifiers-2.13.tar.gz
tar -xzf Class-Method-Modifiers-2.13.tar.gz
cd Class-Method-Modifiers-2.13
perl Makefile.PL PREFIX=$DIR
make
make install
cd ../

wget https://cpan.metacpan.org/authors/id/H/HA/HAARG/Role-Tiny-2.002004.tar.gz
tar -xzf Role-Tiny-2.002004.tar.gz
cd Role-Tiny-2.002004
perl Makefile.PL PREFIX=$DIR
make
make install
cd ../

wget https://cpan.metacpan.org/authors/id/H/HA/HAARG/Sub-Quote-2.006006.tar.gz
tar -xzf Sub-Quote-2.006006.tar.gz
cd Sub-Quote-2.006006
perl Makefile.PL PREFIX=$DIR
make
make install
cd ../

wget https://cpan.metacpan.org/authors/id/H/HA/HAARG/Moo-2.005004.tar.gz
tar -xzf Moo-2.005004.tar.gz
cd Moo-2.005004
perl Makefile.PL PREFIX=$DIR
make
make install
cd ../

wget https://cpan.metacpan.org/authors/id/Y/YA/YANICK/Parallel-ForkManager-2.02.tar.gz
tar -xzf Parallel-ForkManager-2.02.tar.gz
cd Parallel-ForkManager-2.02
perl Makefile.PL PREFIX=$DIR
make
make install
cd ../

echo $LDIR > $DIR/.$NAME.rc

rm -rf Class-Method-Modifiers-2.13 Moo-2.005004 Parallel-ForkManager-2.02 Role-Tiny-2.002004 Sub-Quote-2.006006 Test-Fatal-0.016 Test-Needs-0.002009 Try-Tiny-0.30  Class-Method-Modifiers-2.13.tar.gz Moo-2.005004.tar.gz Parallel-ForkManager-2.02.tar.gz Role-Tiny-2.002004.tar.gz Sub-Quote-2.006006.tar.gz Test-Fatal-0.016.tar.gz  Test-Needs-0.002009.tar.gz Try-Tiny-0.30.tar.gz
