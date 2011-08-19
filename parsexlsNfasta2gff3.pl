#!/homes/dsth/dev/localperl_multithread/bin/perl
use strict;
use Spreadsheet::ParseExcel;
use File::Temp;
use DBI qw(:sql_types);
use Getopt::Long;
##### actually should use cdna2genome - can also use annotation file to give some extra bits?!?

#={{{ alignment parameters
my $exonerate = '/usr/bin/exonerate';
my $gmap = '/homes/dsth/dev/NewCap/gmap/src/gmap';
my $maxintron = 1000000;

my $exonerate_args = ' --model est2genome '
  .'--minintron 20 '
    #.'--query sample-gene.fa --target \
    #aaegypti.SUPERCONTIGS-Liverpool.AaegL1.fa \
  .qq{--maxintron 1000000 $maxintron}
  .'--percent 90 '
  .'--score 100 '
  .'--showvulgar yes '
  .'--softmaskquery no '
  .'--softmasktarget yes '
  .'--showalignment no '
  .'--showtargetgff yes '
  .'--geneseed 250 '
  .'--ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n"';

my $gmap_args = qq{ -K $maxintron }
  .qq{-d Phumanus -D gmapdb/Phumanus }
  #.qq{-d aegypti -D /home/dsth/aegypti }
  .'-f 2';
#=}}}

my $check = 0;
my $dbname = 'dsth_testing_newcap';
my $file = q{};
my $format = q{};
my $sbmid = q{};
my $uploaddir = q{};#'/homes/dsth/dev/NewCap/uploads'.q{/};
my $tmpdir = q{};#'/homes/dsth/dev/NewCap/tempdir/';
my $dbuser = q{};#$ENV{user};
my $dbhost = q{};#$ENV{srv}; 
my $dbpass = q{};#$ENV{pass};
my $dbport = q{};#$ENV{port};

##### ':' makes optional, '=' makes required...
&GetOptions( 
    'check'         =>  \$check,
    'dbname=s'      =>  \$dbname,
    'dbhost=s'      =>  \$dbhost,
    'dbuser=s'      =>  \$dbuser,
    'dbport=s'      =>  \$dbport,
    'dbpass=s'      =>  \$dbpass,
    'file=s'        =>  \$file,
    'format=s'      =>  \$format,
    'sbmid=i'       =>  \$sbmid,
    'tmp_dir=s'     =>  \$tmpdir,
    'upload_dir=s'  =>  \$uploaddir,
);

if ($check) {
    exit(0);
} elsif (!$file     || ($format ne 'fasta' && $format ne 'xls') || !$sbmid || 
         !$tmpdir   || !$uploaddir  || !$dbname     || !$dbhost || !$dbuser||
         !$dbpass   || !$dbport) {
    print qq{\nUsage: $0 -file <filename> -format <xls/fasta> -sbmid <sbmid> 
    -dbhost <hostname> -dbname <db> -dbuser <user> -dbport <port> -dbpass <pass>
    -upload_dir <dir> -tmp_dir <dir>\n};
    exit(0);
}

$file = $uploaddir.'/'.$file;
-e $file || die qq{\ninput file $file does not exist};

#my $tmpgff = File::Temp->new( UNLINK => 0, SUFFIX => '.fa', DIR => '/homes/dsth/dev/NewCap/tempdir' );
my $file_stem = $file;
$file_stem =~ s/.*?([^\/]+)$/$1/;
#print qq{\nusing$file_stem\n};
my $tmp_file = $tmpdir.'/'.$file_stem.'_gff.tmp';
open (my $tmpgff, '>', $tmp_file) || die qq{\ncouldn't write to temp file $tmp_file};
print qq{\nwriting to: $tmp_file};

my @genes;

#={{{ xls handling
if ($format eq 'xls') {

    my $parser   = Spreadsheet::ParseExcel->new();
    my $workbook = $parser->parse($file);

    if ( !defined $workbook ) { die $parser->error(), ".\n"; }

    my $worksheet_n = 0;
    for my $worksheet ( $workbook->worksheets() ) {

        $worksheet_n++;
        #y for now don't bother with sheet 2 - i.e. just template.csv worksheet?!?
        last if ($worksheet_n > 1);

        #my ( $row_min, $row_max ) = $worksheet->row_range();
        #my ( $col_min, $col_max ) = $worksheet->col_range();

        #### do we want to force use of the template - i.e. with the first two rows as title 
        die qq{\nthis does not look like it's based upon the VB template file?!} 
            if (
            $worksheet->get_cell(1,0)->value() ne 'Gene Symbol'
            || $worksheet->get_cell(1,1)->value() ne 'Gene Description'
            || $worksheet->get_cell(1,2)->value() ne 'mRNA Sequence'
            || $worksheet->get_cell(1,3)->value() ne 'Translation Start'
            || $worksheet->get_cell(1,4)->value() ne 'Translation End');

        my $row = 2;

        #y grab entries...
        my ($gs,$gd,$mrna,$ts,$te);

        do {
            $gs = $worksheet->get_cell($row,0)->value();
            $gd = $worksheet->get_cell($row,1)->value();
            $mrna = $worksheet->get_cell($row,2)->value();
            $ts = $worksheet->get_cell($row,3)->value();
            $te = $worksheet->get_cell($row,4)->value();
            $row++;

            #print qq{\n$gs\n$gd\n$mrna\n$ts\n$te\n};

            if (
            $mrna =~ /^[aAcCgGtT]+$/
            && $ts =~ /\d+/ 
            && $te =~ /\d+/
            ) { 
        #        print qq{\nwe have a gene entry: $gs / $gd}; 
        #        print qq{\nhas length }.length($mrna).qq{ with translation start/end = $ts/$te\n};
        #        print qq{>$gs [$gd]\n$mrna\n};
                push @genes, +{ seq => $mrna, 
                    desc => $gd, symbol => $gs, 
                    trnsl_start => $ts, trnsl_end => $te 
                };
            }


        } while ($gs && $gd && $mrna && $ts && $te);
    }

#=}}}    
} else {
#={{{ fasta handling

    my $c = 0;
    open my $ffh, '<', $file || exit(1);
    local $/ = qq{\n>}; 
    while(<$ffh>) {
        if (!$c) {
            s/^>(.*)/$1/ #/ first line will be only one WITH >?!?
        }    
        s/(.*)>$/$1/; #/ all but the last line terminate with \n>
        #print qq{\n*** entry: '$_'\n};
        $c++;
        my ($header,$seq);
        if (/^([^\n]+)\n(.*)/s) {
            $header = $1;
            $seq = $2;
        } else { 
            die qq{\ncould not find header and sequence};
        }
        $seq =~ s/\n//g;

        #### only processing dna?!?
        if ($seq !~ /^[aAcCtTgGnN]+$/) { 
            print qq{\nNOT DNA FASTA!?!};
            exit(1);
        }
        #print qq{\nheader: $header\nseq: $seq};

        my ($symbol,$ts,$te,$desc);
        if ($header =~ /([^;]+);(\d+);(\d+);([^;]+)/) {
            push @genes, +{ desc => $4, symbol => $1, 
              trnsl_start => $2, trnsl_end => $3, seq => $seq, };
        } else {
            push @genes, +{ symbol => $header, seq => $seq };
        }

    }

#=}}}
}

for my $gene (@genes) {
#={{{ process genes one by one with gmap atm...

    #require File::Temp;
    #use File::Temp qw/ :seekable /;

    #$fh = File::Temp->new();
    #$fname = $fh->filename;

    #$fh = File::Temp->new(TEMPLATE => $template);
    #  $fname = $fh->filename;

    my $tmpfasta = File::Temp->new( UNLINK => 1, SUFFIX => '.fa', DIR => '/homes/dsth/dev/NewCap/tempdir' );

    #y overloaded handle?!?
    print $tmpfasta '>'.$gene->{symbol}.qq{\n}.$gene->{seq}.qq{\n};
    #print "Filename is $tmp\n";

    #  $tmp->seek( 0, SEEK_END );

    my @features = &cdna2gff($tmpfasta,$gene,'gmap');

    #/ could load directly or generate gff and then load as per - thus remove those that 
    #/ are identical to current genes too?!?

#=}}}
}

#print qq{\ntemp file: $tmpgff};
#use Log::Log4perl qw(get_logger :levels);
my $dbstr = qq{DBI:mysql:database=$dbname;host=$dbhost;port=$dbport;}
    .qq{user=$dbuser;password=$dbpass};
my $dbhashref = {PrintError => 1, RaiseError => 1};
my $dbi = DBI->connect($dbstr, $dbhashref) or die;
my $sth = $dbi->prepare('update files set tmp_gff = "'.$tmpgff.'" where sbm_id = '.$sbmid);
$sth->execute() || die;

#y should prolly find the start and end of translation sites in gff and confirm its the same as in xls?!?
sub cdna2gff {

    my ($fasta, $gene, $mode) = @_;
    #print qq{$gmap $gmap_args $fasta\n};
    #my @gfflines = `$gmap $gmap_args $fasta`;
    my @gfflines = `$gmap $gmap_args $fasta 2>/dev/null`;
    @gfflines = grep { $_ =~ /^.*([^\t]+\t){8}[^\t]+$/ } @gfflines;
    
    #for my $l (@gfflines) {
    #    $l =~ s/\n//g;
    #    print qq{|\nddd$l}; }
    my $symbol = defined $gene->{symbol} ? ';Symbol="'.$gene->{symbol} : q{"};
    my $desc = defined $gene->{desc} ? ';Desc="'.$gene->{desc} : q{"};
    @gfflines = map { $_ =~ s/\n//; $_.$symbol.$desc.qq{\n} }  @gfflines;
    #print @gfflines;

    #open my $tmpgff, '>', $outfile || die;
    print $tmpgff @gfflines;

    #exit( qq{\nthere's a problem with this gene?!?} if (!@gfflines);
    return 1;
    #unlink 'temp.faI';wq
}

#### factor out the insert routine into a module for use in various procedures?!?

#/ gmap takes cdna sequence and gives traditional gff - i.e. with exon/cds - thus utr etc.,
#/ exonerate gives much more detailed output with just exons - thus would presumably have
#/ to get translation start/stop positions - i.e. cds by translating it?!? - see below

#/ we can of course use the start and end of translations to generate cds?!?

#/ in case of gmap prolly ought to check these values back against those entered?!?
#/ could be uber-lazy and run it twice - i.e. match out start and end first time...

__END__


 gunzip aaegypti.SUPERCONTIGS-Liverpool.AaegL1.fa.gz
 ./gmap-2011-03-28.v3/bin/gmap_setup -d aegypti aaegypti.SUPERCONTIGS-Liverpool.AaegL1.fa

edit makefile locations for fa_coors (utils/bin), gmapindex (src)

 make -f Makefile.aegypti coords
 make -f Makefile.aegypti gmapdb
 make -f Makefile.aegypti install

to output in gff format:

gmap defaults to max intron length of 1000000 (set with -K)
 ./gmap-2011-03-28.v3/src/gmap -d aegypti -D /home/dsth/aegypti sample-gene.fa -f 2

get same result if you set max intron length to same value?!?

exonerate --model est2genome --minintron 20 --query sample-gene.fa --target \
aaegypti.SUPERCONTIGS-Liverpool.AaegL1.fa \
 --maxintron 1000000 \
--percent 90 --score 100 --showvulgar yes --softmaskquery no --softmasktarget yes --showalignment no --showtargetgff yes --geneseed 250 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" > exonerate.gff


# --maxintron 20000 \


>AE [SLC4-like chloride,bicarbonate exchanger]
AGCCAAAAAGTCCATATCTTCTACTACCGTTCGCCGCAAGCTTCCTGCGGGCTTACTAGTGAACGAGATAACAAAGCTAAGTTAAATTTTCACGGTTTTGAGTTGTTCAGAAAAAAAAACAAAGAGGTGTGAGTGGTGATTTGGAGGAAAATTTAGTGTGTGATTTATTTTGTGCAAGAAAAGTTTAATGCTGTTTGAAGAAAAATGAGTCGCCGAGATAATAACGTCAGGAAGCTGTCATTTCTTGGTTTCAACACTAAAGAGACAAGCAACGATGATCCCAACCACGTGCAGTTGGACGATGAGATGGAAAAGGTGTTCGGATCGGTTGGGACCGATAAGGAGCGGTTCGAACTGAAGCGCCTCAACGATGAGGTGCTGGTGGACAACTCGCCGTTAAAGTACGACGAATCACACCGATCCGTGGATAATCGTCCGTTACTGTCCTCGGCGGCCCTTCGACCGATGGATAGCCCCGCTGCAGGCGGCGGAGGAGCAGGAGGTGGAGGAGGTGGTGGGCCAAGCAGTGGACAGCAGAACTCTCCACCGACGAGCACCCCAACGAGTCCTATTTCGATTTCATCCAAAAGCGAAACTACCACCAGACCATCGCATGATACGACGCTGGCGGACAACACATCCAACGATTTCAGCGAAGCTGTCCAAGATGAACCTGTGGTGGACCAGGGAACGCTGCAAGGGGAACAGTGGGATGCGAGCGCCAGGAGAAACGTTCATTTCGACAGCAAGGACCGTCCCCCACAGTTCGAGGGGCTTCAGATTGAAGACAGCAATGAGGAACGTCGACGCCGTACCGAGCGTCATCATCCTCATAAATCGAGGAAGTTCTCCTTGCAAGAGTACCATCCTGAGTGGAGGCGTCAAAGTGGAGCAGAAGGCGCCAGCACCACTACCAGGAGGGTTTCCGTGCAACCCGAGGATGCCACATTACAGGAAGCCGACATTGACGAATTGACCTCTCACCGGTCGGATGACCCTCGAGCGTTGCGCCGTCACAAGGTAAGCGCCCAGTCGCAAACGCAGCCGGGCGGTCCTTCGATGGTCAACATCAACCGGAAGGATGGTGACAAACTGCAACATCTGCTTCCGTCGAACAAATTCAAGAAAATGTATGACCACAGCCCCCATGAGGTGTTCGTTCAGCTGGACGAACTTACCGGATCCGGCGAAGACCGGGAATGGAAGGAAACGGCCCGTTGGATCAAGTATGAGGAAGATGTGGAAGAAGGCGCCGATCGCTGGGGTCGTCCACATGTGGCTTCGCTGTCGTTCCATTCGCTGCTAAATCTTCGTCGTTGCCTCGAGACCGGAGTGGTGTTGATGGACTTGGAGGAAAAGGATCTGCCCTCCGTTGCTTACCGGATTGTGGAACAGATGGTCATCGACGAGCTGATCCATGAGGATGACAAACCTACGATTATGCGAGCTTTGTTGCTCAGACATCGTCACGTGAACGAACACTCGCACGGTGGATTCCATTTCGGACCAAAGCGGAAGTACAGCAGCTACAGTAGCTTACAGAGCGTGGACGACAAGAAGCCACGAATTGTGCCATCGAGTGAGATTAACGGCCACGGACATGGCGAGACCAAGATCAACATGCACGAAGAAACGTACACCTCGTCACAGGAGGACATCAAAATGCGCACCCAGAAAGAATCGATCCTGAAGCGTATACCCGAAGGTGCCGAAGCCACAACAGTGCTGGTAGGATCGGTTGACTTTCTGGAGCAACCAACGATCGCTTTTGTCCGCTTGGCAGAGGGTATTCCGATGCCTAGTATTACGGAGGTGCCAATTCCGGTACGCTTCCTGTTCATTCTGTTGGGACCAAAGACTGCCGAGCTGGACTACCATGAAGTGGGTCGATCGATTGCCACACTTATGTCCAACGAGCATTTCCATGATATTGCATACAGGGCAGACGATAGGAAGGACTTACTGTCGGCTATTAATGAATTTTTGGACGATTCGATTGTGCTGCCACCTGGCAAGTGGGAACGGCAGGCTTTGTTGCCGTTCGACGAGTTGAAGGCCAAGAGTGATATGATCCGACTGAGAAAGAAGAAAGCAATCGATGAAAAGATCAAGAGCAAACAGCCACTACTGACCAGTGAAGAGGAGAAGAAGCTTTTGGCGGCAGCTGAAGGAGATGGGAAAAAGCCTACTAAAAATCCTTTGGAGAAGACCCACCGATTATGGGGTGGTCTCATAAACGATATTAAAAGGAGGTATCCGATGTATAAGAGCGATATTAAAGATGGATTAAACACGGAAACTCTTGCAGCTACTGTATTCATGTACTTTGCTGCCTTATCAACGGCCATTACATTTGGAGGTCTTGCCTCAGATAAAACACACAATTTGATTGGAATTTCTGAAACTTTGGTATCTGCGTCGATGGTTGGCGTTGTGTTTCATTTGTTTTCCGGTCAACCGTTGGTCATCATCGGTACCACTGGTCCTCTATTGCTATTCGATGAAGCGCTCAATCAATTTTGTATCTCAAACAATTTCAGTTTCTTGACAGTGCGAGTCTACGTAGGTATTTGGCTAGCCGTTATTGCTTTAGTTGTGTCTGCCTTTGAAGGTAGTGTCTACGTACGGCTATTCACTCGATTCACTCAAGAAATCTTTTCCGCTTTGATCACATTGATCTATATCGTCGAAACTGTTATGAAGCTGGTCTACGTTTATGGCAGACATCCTCTTTTGGCAGAGTATCAATACAAAAACTTGACCATCACACCGCAGCCTTTGCTCGAAAGCGAATCCAACTCAACAATTGATTCTTCGGTTGGCCTCCTAGCGGAAAGCCTTACGGCAGTTGCAAATGCTACCGTAGCACCGTTATCTGATAACTTACTGCCTCCATCGGATGTCGTTGGACCATTGAATCAACCAAACACAGCCTTATTCTGCACAATTCTAACCTTGGGAACATTCTCGCTGGCTTACTATCTGAAATTGTTCAGAAACTCTCACTTCCTGGGACGTAACGCTCGACGTGCCCTTGGTGATTTCGGCGTGCCAATCTCGATTGCTATCTTCGTGCTCATCGATTATATGATTCCGCAAGTATACACCGAAAAGCTCAGCGTTCCGGAGGGTCTTTCTCCGAGTGACGAAACCCGTCGAGGATGGATCATTCCTTTGGGTGGAGTCCCTAGCTGGATGCCATTCATGGCCGGTATCCCTGCCCTGTTGGTATACATCTTGATCTTCATGGAGACCCACATTTCCGAGCTGATCGTAGACAAACCGGAACGCGGACTAAAGAAAGGATCGGGTCTGCACATGGACATCGTCCTGCTGTGCTTCCTCAATACGGTGTGTGGATTCTTTGGAATGCCGTGGCATTGTGCTGCTACCGTACGATCGGTAACACACGTATCCGCAGTGACGATCATGTCTAGAACTCATGCACCTGGTGAATCTCCACACATTACGGACGTTAAGGAGCAGCGTATCTCCGGATTCTTCGTTTCGCTGATGGTGGGCCTTTCGGTCACTATGGCACCCATTCTGCGACTCATCCCTATGTCCGTGCTGTTTGGCGTTTTCCTGTACATGGGTATCGCTTCCATGAGTGGAGTTCAATTCTTCGAAAGATTGCGACTATTCTTCATGCCCGTAAAGCATCACCCGCAGGTGCCGTTTGTGCGTCGGGTGCCTTCCTGGAAGATGCACATCTTTACGGCGGTGCAGATTCTTGCCTTGGCCATGTTGTGGGCCGTCAAGTCATCGGCCTTCTCTTTGGCGTTCCCGTTCTTCCTCATCATGATGGTGCCAATCAGGAAGCAGATGGAAAGAATTTTCAGCCCATTGGAACTCCGAGCGCTGGACGGCAGCCAACCCAACGAAGGAGCCGAAGATGAACCTGACTTCTACGAACAAGCGCCAATTCCTGCTTAAGAACTTAAGCAGCGTAGTATTAGGTGAAAAGTGCAACCAACATAAACCAATCCAATCTGTATTGTGTATTGATTTTGCGACGTTTATTATGTTTGAAAAAGCCGGAACACGTTATTCGTCACGATTGTTTCACTTATTTACAACTTGTATGGTAGAATTAGTGCAAAATTAAGAATGCAAAATCAATTCGTCCATGCATATTATACATTGTTTTCAGAACAGCTCGATCGGACGAAACTTCTATGAAGTACAATAAATTTTTGTTAAATGTTTACCAGTTACTCATCGTGCACACAAATTTAGTGAGGTTGAGCGGTTAGAGTCTATACTTATATTTTTCTACGAATGTGTTTTATTAGTCACATGACTAGGATAACACCATTGAATAGAAACAATAAATTTCAAATAAAATTGATACTTTCTCTACGCTAAATTTTCCCTTTGAATGAAATTCATACAAAAAAAAAAAAAAAAAAAAA


supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	gene	303213	439809	.	+	.	ID=AE.path1;Name=AE
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	mRNA	303213	439809	.	+	.	ID=AE.mrna1;Name=AE;Parent=AE.path1;Coverage=99.8;Identity=99.6
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	303213	303295	98	+	.	ID=AE.mrna1.exon1;Name=AE;Parent=AE.mrna1;Target=AE 1 83 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	303739	303919	100	+	.	ID=AE.mrna1.exon2;Name=AE;Parent=AE.mrna1;Target=AE 84 264 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	416149	416595	100	+	.	ID=AE.mrna1.exon3;Name=AE;Parent=AE.mrna1;Target=AE 265 711 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	424776	424882	100	+	.	ID=AE.mrna1.exon4;Name=AE;Parent=AE.mrna1;Target=AE 712 818 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	425351	425486	100	+	.	ID=AE.mrna1.exon5;Name=AE;Parent=AE.mrna1;Target=AE 819 954 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	425550	425990	100	+	.	ID=AE.mrna1.exon6;Name=AE;Parent=AE.mrna1;Target=AE 955 1395 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	426053	426199	100	+	.	ID=AE.mrna1.exon7;Name=AE;Parent=AE.mrna1;Target=AE 1396 1542 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	436473	438367	99	+	.	ID=AE.mrna1.exon8;Name=AE;Parent=AE.mrna1;Target=AE 1543 3437 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	438430	438624	100	+	.	ID=AE.mrna1.exon9;Name=AE;Parent=AE.mrna1;Target=AE 3438 3632 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	438695	438923	98	+	.	ID=AE.mrna1.exon10;Name=AE;Parent=AE.mrna1;Target=AE 3633 3861 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	exon	439264	439809	99	+	.	ID=AE.mrna1.exon11;Name=AE;Parent=AE.mrna1;Target=AE 3862 4407 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	CDS	303860	303919	100	+	0	ID=AE.mrna1.cds1;Name=AE;Parent=AE.mrna1;Target=AE 205 264 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	CDS	416149	416595	100	+	0	ID=AE.mrna1.cds2;Name=AE;Parent=AE.mrna1;Target=AE 265 711 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	CDS	424776	424882	100	+	0	ID=AE.mrna1.cds3;Name=AE;Parent=AE.mrna1;Target=AE 712 818 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	CDS	425351	425486	100	+	2	ID=AE.mrna1.cds4;Name=AE;Parent=AE.mrna1;Target=AE 819 954 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	CDS	425550	425990	100	+	0	ID=AE.mrna1.cds5;Name=AE;Parent=AE.mrna1;Target=AE 955 1395 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	CDS	426053	426199	100	+	0	ID=AE.mrna1.cds6;Name=AE;Parent=AE.mrna1;Target=AE 1396 1542 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	CDS	436473	438367	99	+	0	ID=AE.mrna1.cds7;Name=AE;Parent=AE.mrna1;Target=AE 1543 3437 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	CDS	438430	438624	100	+	2	ID=AE.mrna1.cds8;Name=AE;Parent=AE.mrna1;Target=AE 3438 3632 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	CDS	438695	438923	98	+	2	ID=AE.mrna1.cds9;Name=AE;Parent=AE.mrna1;Target=AE 3633 3861 +
supercontig_AaegL1_supercont1.28_1_3768427_1	aegypti	CDS	439264	439336	100	+	0	ID=AE.mrna1.cds10;Name=AE;Parent=AE.mrna1;Target=AE 3862 3934 +
###


Command line: [exonerate --model est2genome --minintron 20 --query sample-gene.fa --target aaegypti.SUPERCONTIGS-Liverpool.AaegL1.fa --maxintron 1000000 --percent 90 --score 100 --showvulgar yes --softmaskquery no --softmasktarget yes --showalignment no --showtargetgff yes --geneseed 250]
Hostname: [dsth-laptop]
vulgar: AE 0 4415 + supercontig:AaegL1:supercont1.28:1:3768427:1 303212 439819 + 21795 M 83 83 5 0 2 I 0 439 3 0 2 M 181 181 5 0 2 I 0 112225 3 0 2 M 447 447 5 0 2 I 0 8176 3 0 2 M 107 107 5 0 2 I 0 464 3 0 2 M 136 136 5 0 2 I 0 59 3 0 2 M 441 441 5 0 2 I 0 58 3 0 2 M 147 147 5 0 2 I 0 10269 3 0 2 M 1895 1895 5 0 2 I 0 58 3 0 2 M 195 195 5 0 2 I 0 66 3 0 2 M 229 229 5 0 2 I 0 336 3 0 2 M 540 540 G 0 2 M 14 14
# --- START OF GFF DUMP ---
#
#
##gff-version 2
##source-version exonerate:est2genome 2.2.0
##date 2011-07-04
##type DNA
#
#
# seqname source feature start end score strand frame attributes
#
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	gene	303213	439819	21795	+	.	gene_id 1 ; sequence AE ; gene_orientation +
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	utr5	303213	303295	.	+	.	
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	303213	303295	.	+	.	insertions 0 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice5	303296	303297	.	+	.	intron_id 1 ; splice_site "GT"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	intron	303296	303738	.	+	.	intron_id 1
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice3	303737	303738	.	+	.	intron_id 0 ; splice_site "ag"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	utr5	303739	303919	.	+	.	
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	303739	303919	.	+	.	insertions 0 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice5	303920	303921	.	+	.	intron_id 2 ; splice_site "GT"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	intron	303920	416148	.	+	.	intron_id 2
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice3	416147	416148	.	+	.	intron_id 1 ; splice_site "AG"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	utr5	416149	416595	.	+	.	
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	416149	416595	.	+	.	insertions 0 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice5	416596	416597	.	+	.	intron_id 3 ; splice_site "GT"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	intron	416596	424775	.	+	.	intron_id 3
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice3	424774	424775	.	+	.	intron_id 2 ; splice_site "AG"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	utr5	424776	424882	.	+	.	
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	424776	424882	.	+	.	insertions 0 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice5	424883	424884	.	+	.	intron_id 4 ; splice_site "GT"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	intron	424883	425350	.	+	.	intron_id 4
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice3	425349	425350	.	+	.	intron_id 3 ; splice_site "AG"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	utr5	425351	425486	.	+	.	
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	425351	425486	.	+	.	insertions 0 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice5	425487	425488	.	+	.	intron_id 5 ; splice_site "GT"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	intron	425487	425549	.	+	.	intron_id 5
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice3	425548	425549	.	+	.	intron_id 4 ; splice_site "AG"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	utr5	425550	425990	.	+	.	
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	425550	425990	.	+	.	insertions 0 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice5	425991	425992	.	+	.	intron_id 6 ; splice_site "GT"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	intron	425991	426052	.	+	.	intron_id 6
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice3	426051	426052	.	+	.	intron_id 5 ; splice_site "AG"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	utr5	426053	426199	.	+	.	
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	426053	426199	.	+	.	insertions 0 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice5	426200	426201	.	+	.	intron_id 7 ; splice_site "GT"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	intron	426200	436472	.	+	.	intron_id 7
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice3	436471	436472	.	+	.	intron_id 6 ; splice_site "AG"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	utr5	436473	438367	.	+	.	
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	436473	438367	.	+	.	insertions 0 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice5	438368	438369	.	+	.	intron_id 8 ; splice_site "GT"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	intron	438368	438429	.	+	.	intron_id 8
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice3	438428	438429	.	+	.	intron_id 7 ; splice_site "AG"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	utr5	438430	438624	.	+	.	
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	438430	438624	.	+	.	insertions 0 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice5	438625	438626	.	+	.	intron_id 9 ; splice_site "GT"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	intron	438625	438694	.	+	.	intron_id 9
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice3	438693	438694	.	+	.	intron_id 8 ; splice_site "AG"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	utr5	438695	438923	.	+	.	
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	438695	438923	.	+	.	insertions 0 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice5	438924	438925	.	+	.	intron_id 10 ; splice_site "GT"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	intron	438924	439263	.	+	.	intron_id 10
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	splice3	439262	439263	.	+	.	intron_id 9 ; splice_site "AG"
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	exon	439264	439819	.	+	.	insertions 2 ; deletions 0
supercontig:AaegL1:supercont1.28:1:3768427:1	exonerate:est2genome	similarity	303213	439819	21795	+	.	alignment_id 1 ; Query AE ; Align 303213 1 83 ; Align 303739 84 181 ; Align 416149 265 447 ; Align 424776 712 107 ; Align 425351 819 136 ; Align 425550 955 441 ; Align 426053 1396 147 ; Align 436473 1543 1895 ; Align 438430 3438 195 ; Align 438695 3633 229 ; Align 439264 3862 540 ; Align 439806 4402 14
# --- END OF GFF DUMP ---
#
-- completed exonerate analysis


