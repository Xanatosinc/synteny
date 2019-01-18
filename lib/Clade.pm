package Clade;

use strict;
use warnings;

sub new {
    my $class = shift;
    my $self = {
        genome => shift,
    };
    bless $self, $class;
    return $self;
}

1;
