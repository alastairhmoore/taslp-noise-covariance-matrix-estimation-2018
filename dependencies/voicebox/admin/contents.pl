# Routine to create contents file for matlab
#read in m-file names in the parent directory

while (<..\\*.m>) {
    tr/A-Z.\\/a-z/d;
    chop;
    if ( !/contents/ && !/(lpc|rot)..2../ ) {
        $index{$_} = 0;
    }
}
$unc = 0;
%mis = 0;
if ( open( CFILE, "..\\Contents.m" ) ) {
    @lines = <CFILE>;
    open( CFILE, ">..\\Contents.m" ) || die "Cannot write to contents.m\n";
    $st = 1;
    foreach (@lines) {
        if ( $st == 1 ) {
            if (/^%%%%%%%%%%%/) {
                print CFILE '%' x 80, "\n";
                $st = 2;
                while ( ( $mfile, $stat ) = each %index ) {
                    if ( $stat == 0 && $mfile !~ /(peak2dquad)/ ) {
                        if ( $unc == 0 ) {
                            $unc = 1;
                            print CFILE "% === Unclassified ===\n";
                        }
                        print CFILE "%   $mfile", ' ' x ( 14 - length($mfile) ),
                          "- \n";
                    }
                }
                while ( ( $mfile, $stat ) = each %index ) {
                    if ( $stat == 2 && $mfile !~ /(lpc|rot)..2../ ) {
                        if ( $mis == 0 ) {
                            $mis = 1;
                            print CFILE "% === Missing ===\n";
                        }
                        print CFILE "%   $mfile\n";
                    }
                }
            }
            else {
                print CFILE $_;
                if (/^%   (\s*)(\w+) /) {
                    if ( length($1) < 10 ) {
                        if ( exists $index{$2} ) {
                            $index{$2} = 1;
                        }
                        else {
                            $index{$2} = 2;
                        }
                    }
                }
            }
        }
        elsif ( $st == 2 ) {
            if (/^%\s*Copyright/) {
                $st = 3;
                print CFILE "\n", '%' x 80, "\n$_";
            }
        }
        else {
            print CFILE $_;
        }
    }
}
else {
    open( CFILE, ">..\\Contents.m" ) || die "Cannot write to contents.m\n";
    print CFILE "% Voicebox: speech processing toolbox\n";
    print CFILE "% Unclassified\n";
    while ( ( $mfile, $stat ) = each %index ) {
        print CFILE "%   $mfile\n";
    }
    print CFILE "% Copyright (c) Mike Brookes 1998\n";
}
close(CFILE);
system "notepad ..\\Contents.m";
