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
    my $usage = "Usage $0 [--debug] [--caller_list] [--caller_delta] --clades CLADES_FILE --clusters GENE_CLUSTER_FILE\n";
    my ( $debug_opt, $caller_list_opt, $caller_delta_opt, $cladeFilename, $clusterFilename );
    GetOptions(
        'debug' => \$debug_opt,
        'caller_list' => \$caller_list_opt,
        'caller_delta' => \$caller_delta_opt,
        'clades=s' => \$cladeFilename,
        'clusters=s' => \$clusterFilename,
    ) or die $usage;
    if ($cladeFilename && $clusterFilename) {
        my %clades = scanClades($cladeFilename, $debug_opt);
        my %clusters = scanClusters(\%clades, $clusterFilename, $debug_opt);

        print " ### FILTERING ###\n" if $debug_opt;
        printf " /--------------------| %-8s | %-8s | %-8s | %-8s |\n", "GENOME", "DELTA", "CALLER1", "CALLER2" if $caller_delta_opt;
        printf "| %-8s | %-8s | %-8s | %-8s | %-8s | %-8s |\n", "CLADE", "GC_COUNT", "COUNT", "MIN", "MAX", "AVG" if ($debug_opt || $caller_delta_opt);
        printf "%s\t%s\t%s\t%s\n", "Clade", "Genome_Name", "Caller_Id", "Cluster_Id" if $caller_list_opt;
        foreach my $cladeId (sort keys %clades) {
            #print "[$cladeId]\n" if $debug_opt;
            my %filteredClusters;

            # Check each GC; if it contains all of the clade's genomes, keep it
            foreach my $clusterCheckId (sort keys %clusters) {
                if ($clusters{$clusterCheckId}->checkGenomes(keys $clades{$cladeId})) {
                   $filteredClusters{$clusterCheckId} = $clusters{$clusterCheckId};
                }
            }
            my $cluster_per_clade_count = scalar keys %filteredClusters;

            # Arrange data such that %genomeCallers{genome_id}{caller_id} = cluster_id
            my %genomeCallers;
            foreach my $genome (keys $clades{$cladeId}) {
                foreach my $clusterId (sort keys %filteredClusters) {
                    my @callers = keys $filteredClusters{$clusterId}->{genomes}{$genome};
                    foreach (@callers) {
                        $genomeCallers{$genome}{$_} = $clusterId;
                    }
                }
            }
            my $calls_per_clade_count = scalar keys %genomeCallers;
            my $min_calls_per_genome_count = 9999;
            my $max_calls_per_genome_count = 0;
            my $sum_calls_per_clade = 0;
            foreach my $genome (sort keys %genomeCallers) {
                my $caller_count = int( scalar keys $genomeCallers{$genome} );
                $max_calls_per_genome_count = max(
                    $caller_count,
                    $max_calls_per_genome_count
                );
                $min_calls_per_genome_count = min($caller_count, $min_calls_per_genome_count);
                $sum_calls_per_clade += $caller_count;

                my ($max_caller_delta, $prev_caller, $next_caller, $max_prev_caller, $max_next_caller) = (0, 0, 0, 0, 0);

                foreach my $callerId (sort keys $genomeCallers{$genome}) {
                    $prev_caller = $next_caller;
                    $next_caller = $callerId;
                    if ($prev_caller) {
                        if ($next_caller - $prev_caller > $max_caller_delta) {
                            $max_caller_delta = $next_caller - $prev_caller;
                            $max_prev_caller = $prev_caller;
                            $max_next_caller = $next_caller;
                        }
                    }

                    my $clusterId = $genomeCallers{$genome}{$callerId};
                    printf "%s\t%s\t%s\t%s\n", $cladeId, $genome, $callerId, $clusterId if $caller_list_opt;
                }
                outputDeltaInfo($genome, $max_caller_delta, $max_prev_caller, $max_next_caller) if $caller_delta_opt;

            }
            my $avg = $sum_calls_per_clade / $calls_per_clade_count;
            printf "| %-8s | %8.0f | %8.0f | %8.0f | %8.0f | %8.2f |\n\n", (
                $cladeId,
                $cluster_per_clade_count,
                $calls_per_clade_count,
                $min_calls_per_genome_count,
                $max_calls_per_genome_count,
                $avg,
            ) if $debug_opt || $caller_delta_opt;
        }
    } else {
        print $usage;
    }
}

sub outputDeltaInfo {
    my ( $genome, $max_caller_delta, $max_prev_caller, $max_next_caller ) = @_;

    printf " /--------------------| %-8s | %-8s | %-8s | %-8s |\n", (
        $genome, $max_caller_delta, $max_prev_caller, $max_next_caller
    );
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
