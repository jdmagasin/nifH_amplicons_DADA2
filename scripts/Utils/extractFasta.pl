#!/usr/bin/perl

#
# Extract from a FASTA the sequences specified in an ID's file (one ID per line).
# ID's are everything between the '>' and the first whitespace (or newline).
# The FASTA can be gzip'd or not.  Ignore comment lines (which begin with ';')
# and blank lines.
#

use strict;
use warnings;

my ($ids, $fasta) = @ARGV;
defined($ids) and -f $ids or die('IDs file is...?');
defined($fasta) and -f $fasta or die('fasta file is...?');

my %wanted;  # Keys are ids of sequences wanted.
my @order;
open(FH,"<$ids") or die();
while (<FH>) {
    next if (/^;/);
    chomp;
    $wanted{$_} = '';  # Set to sequence when found
    push(@order,$_);
}
close(FH);


# Get sequences wanted. Extremely simple, no error checking (in FASTAs
# or for duplicates) ID is considered everything up to first WS of a
# defline.
my $id = undef;
my $cmd = ($fasta =~ /\.gz$/ ? "cat $fasta | gunzip |" : "<$fasta");
my %deflines;
open(FH, $cmd) or die();
while (<FH>) {
    next if (/^\s*$/);
    chomp;
    if (/>([^\s]+)/) {  # /w+ doesn't get periods (at least)
        $id = (defined($wanted{$1}) ? $1 : undef);
        defined($id) and ($deflines{$id} = $_);
    } elsif (defined($id)) {
        $wanted{$id} = $wanted{$id}.$_;
    }
}
close(FH);

foreach (@order) {
    print STDOUT "$deflines{$_}\n$wanted{$_}\n\n" if (defined($deflines{$_}));
}
foreach (@order) {
    print STDERR "Could not find $_\n" if (!defined($wanted{$_}));
}
