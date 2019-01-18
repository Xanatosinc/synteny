#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long qw(GetOptions);
use List::MoreUtils qw(any);
use List::Util qw(min max);
use Number::Format qw(format_number);

use lib::GeneCluster;

sub main {
    my $usage = "Usage $0 --debug --caller_list --clades CLADES_FILE --clusters GENE_CLUSTER_FILE\n";
    my ( $debug_opt, $caller_list_opt, $cladeFilename, $clusterFilename );
    GetOptions(
        'debug' => \$debug_opt,
        'caller_list' => \$caller_list_opt,
        'clades=s' => \$cladeFilename,
        'clusters=s' => \$clusterFilename,
    ) or die $usage;
    if ($cladeFilename && $clusterFilename) {
        my %clades = scanClades($cladeFilename, $debug_opt);
        my %clusters = scanClusters(\%clades, $clusterFilename, $debug_opt);

        print " ### FILTERING ###\n" if $debug_opt;
        printf "| %-8s | %-8s | %-8s | %-8s | %-8s | %-8s |\n", "CLADE", "GC_COUNT", "COUNT", "MIN", "MAX", "AVG" if $debug_opt;
        printf "%s\t%s\t%s\n", "Clade", "Genome_Name", "Caller_Id" if $caller_list_opt;
        foreach my $cladeId (sort keys %clades) {
            #print "[$cladeId]\n" if $debug_opt;
            my %filteredClusters;

            # Check each GC; if it contains all of the clade's genomes, keep it
            foreach my $clusterCheckId (sort keys %clusters) {
                if ($clusters{$clusterCheckId}->checkGenomes(keys $clades{$cladeId})) {
                   $filteredClusters{$clusterCheckId} = $clusters{$clusterCheckId};
                }
            }
            my $filterCount = scalar keys %filteredClusters;

            my %genomeCallers;
            foreach my $genome (keys $clades{$cladeId}) {
                foreach my $clusterId (sort keys %filteredClusters) {
                    my @callers = keys $filteredClusters{$clusterId}->{genomes}{$genome};
                    foreach (@callers) {
                        $genomeCallers{$genome}{$_} = $clusterId;
                    }
                }
            }
            my $count = scalar keys %genomeCallers;
            my $min = 9999;
            my $max = 0;
            my $sum = 0;
            foreach my $genome (sort keys %genomeCallers) {
                my $callerCount = int( scalar keys $genomeCallers{$genome} );
                $max = max(
                    $callerCount,
                    $max
                );
                $min = min($callerCount, $min);
                $sum += $callerCount;

                #printf "%s callerIds\n", $callerCount if $debug_opt;
                foreach my $callerId (sort keys $genomeCallers{$genome}) {
                    #printf "%-8s : %s\n", $callerId, $genomeCallers{$genome}{$callerId} if $debug_opt;
                    printf "%s\t%s\t%s\n", $cladeId, $genome, $callerId if $caller_list_opt;
                }
            }
            my $avg = $sum / $count;
            printf "| %-8s | %8.0f | %8.0f | %8.0f | %8.0f | %8.2f |\n", $cladeId, $filterCount, $count, $min, $max, $avg if $debug_opt;
        }
    } else {
        print $usage;
    }
}

sub scanClades {
    my $filename = shift;
    my $debug = shift;

    my %clades;
    my %uniqueGenomes;
    
    open(my $cladeFile, '<', $filename)
        or die "Could not open '$filename' $!";
    
    my $header = <$cladeFile>;
    while (<$cladeFile>) {
        chomp;
        my @row = split('\t');
        $clades{$row[1]}{$row[0]} = 1;
        $uniqueGenomes{$row[0]} = 1;
    }
    if ($debug) {
        print " ##CLADES##  ##GENOMES##\n";
        for (sort keys %clades) {
            printf "[ %-8s ] %s\n", $_, join(", ", sort keys $clades{$_});
        }
        printf "%s unique genomes scanned from clade file\n", format_number(scalar keys %uniqueGenomes);
        printf "[ %s ]\n", join(", ", sort keys %uniqueGenomes);
        print "\n";
    }
    return %clades;
}

sub scanClusters {
    my %clades = %{$_[0]};
    my $filename = $_[1];
    my $debug = $_[2];

    my %clusters;
    my %uniqueGenomes;

    open(my $clusterFile, '<', $filename)
        or die "Could not open '$filename' $!";

    my $header = <$clusterFile>;
    while (my $line = <$clusterFile>) {
        chomp $line;
        my @row = split('\t', $line);
        my ( $clusterId, $genomeName, $callerId ) = ( $row[1], $row[3], $row[4] );
        $uniqueGenomes{$genomeName} = 1;
        if ($clusters{$clusterId}) {
            $clusters{$clusterId}->{genomes}{$genomeName}{$callerId} = 1;
        } else {
            my %genomes = (
                $genomeName => {
                    $callerId => 1,
                },
            );
            $clusters{$clusterId} = GeneCluster->new($clusterId, \%genomes);
        }
    }
    if ($debug) {
        print " ### CLUSTERS ###\n";
        printf "%s clusters scanned.\n", format_number(scalar keys %clusters);
        my $firstClusterId = (sort keys %clusters)[0];
        my $firstCluster = $clusters{$firstClusterId};
        print "\tFirst cluster ";
        $firstCluster->summarize();
        my $firstGenomeId = (sort keys $firstCluster->{genomes})[0];
        printf "\t\tFirst Genome: %s # Caller IDs: %s\n", $firstGenomeId, scalar keys $firstCluster->{genomes}{$firstGenomeId};
        my $firstCallerId = (sort keys $firstCluster->{genomes}{$firstGenomeId})[0];
        printf "\t\t\tFirst Caller: %s\n", $firstCallerId;
        printf "%s unique genomes scanned from cluster file\n", format_number(scalar keys %uniqueGenomes);
        printf "[ %s ]\n", join(", ", sort keys %uniqueGenomes);
        print "\n";
    }

    return %clusters;
}

&main();
