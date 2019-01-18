package GeneCluster;

use strict;
use warnings;

use Data::Dumper;

sub new {
    my $class = shift;
    my $self = {
        cluster_id => shift,
        genomes => shift,
    };
    bless $self, $class;
    return $self;
}

sub summarize {
    my ( $self ) = @_;
    printf "Id: %-12s # Genomes: %s\n", $self->{cluster_id}, scalar keys $self->{genomes};
}

sub checkGenomes {
    my ( $self, @genomes ) = @_;
    foreach my $gene (@genomes) {
        if ( !exists( $self->{genomes}{$gene} ) ) {
            return 0;
        }
    }
    return 1;
}

1;
