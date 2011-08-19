#!/homes/dsth/dev/localperl_multithread/bin/perl
use strict;
use warnings;
use IO::Handle;
use File::Map 'map_file';
use Log::Log4perl qw(get_logger :levels);
use threads ('yield',
 'stack_size' => 64*4096,
 'exit' => 'threads_only',
 'stringify', 
);
use threads::shared;
use DBI qw(:sql_types);
use Getopt::Long;

my $modified :shared;
my $features :shared; # my $var :shared; # my %hsh :shared; # my @ary :shared;
$modified = 0;
$features = 0;

#y log4 set-up
my $log4 = get_logger('Cap::Gff');
$log4->level($DEBUG);
my $appender1 = Log::Log4perl::Appender->new('Log::Dispatch::Screen');
my $layout = Log::Log4perl::Layout::PatternLayout->new("CAPPROC : [%p] : %d : %m%n");
#my $layout = Log::Log4perl::Layout::PatternLayout->new("[%p] CAPPROC : %d $0 > %m%n");
#my $layout = Log::Log4perl::Layout::PatternLayout->new("%p %d %M > %m%n");
$appender1->layout($layout);
$log4->add_appender($appender1);

#={{{ a little debug
my $crc_debug = q{};
my $cd = q{};
my $cd1 = q{};
my $cd2 = q{};
my $cd3 = q{};
my $cd4 = q{};
if ($crc_debug) { ($cd, $cd1, $cd2, $cd3, $cd4) = qw/_=_ 1@ 2@ 3@ 4@/; }
#=}}}

my $schema = 0;
my $interval = 35;
my $db = q{};
my $dbuser = q{};#$ENV{user};
my $dbhost = q{};#$ENV{srv}; 
my $dbpass = q{};#$ENV{pass};
my $dbport = q{};#$ENV{port};
my $crc64 = '/homes/dsth/dev/crc64';
my $table = 'files';
my $threads = 17;
my $upload_dir = q{};#q{/homes/dsth/dev/NewCap/uploads};
my $name_length = 256;
my $prohibited_types = qr{three_prime_UTR|five_prime_UTR|protein|polypeptide|contig|supercontig};
my $insert_only = 0; #/ dev
my ($user,$sbmid,$uid,$sbm_file);
my $initial_populate = 0;
my @permitted_biotypes = qw{mRNA miRNA pseudogene ncRNA};
my $mode = q{};
my $file = q{};
my $check = 0;

=head1 SYNOPSIS

# load schema:
# perl gffCapProcessing.pl schema | mysql -udruser -pdruser cap

# extract gff - i.e. get rid of fasta entries
# perl -0ne 'my @gff;s/(##gff.+?)##FASTA/push @gff, $1/esmg; print @gff' with_fasta.gff > gff_only.gff;

# remove single features - i.e. nothing to do with genes... - but break rule of genes being only features with ids and without parents?!?
# grep -v direct_repeat ~/dev/Imports/WTSIscripted/Smansoni/merged.gff > Smansoni_woDirectRepeat.gff
# grep -v nucleotide_match  Smansoni_woDirectRepeat.gff > t; mv t Smansoni_woDirectRepeat.gff

# set-up-crc64-table: 
# perl gffCapProcessing.pl initialise <filename>

# run-an-import: 
# perl gffCapProcessing.pl

# recover gff:
# mysql -udruser -pdruser cap -B -e 'select seqid, source, type, start, end, score, strand, phase, annot from feat'

# recover name->annot for xref parsing?
# mysql -udruser -pdruser cap -B -e 'select name, annot from feat where name is not null'

=cut

##### ':' makes optional, '=' makes required...
&GetOptions( 
    'check'         =>  \$check,
    'schema'        =>  \$schema,
    'initialise'    =>  \$initial_populate,
    'crc64=s'       =>  \$crc64,
    'mode=s'        =>  \$mode,
    'file=s'        =>  \$sbm_file,
    'dbhost=s'      =>  \$dbhost,
    'dbname=s'      =>  \$db,
    'dbuser=s'      =>  \$dbuser,
    'dbport=s'      =>  \$dbport,
    'dbpass=s'      =>  \$dbpass,
    'upload_dir=s'  =>  \$upload_dir,
    'capuser=s'     =>  \$user,
    'capuid=s'      =>  \$uid,
    'sbmid=s'       =>  \$sbmid,
);

if ($check) {
    exit(0);
} elsif ($schema) {
    &schema();
    exit(0);
} elsif ($initial_populate) {
    if (!$sbm_file || !$db || !$dbhost || !$dbuser|| !$dbpass || !$dbport) {
        print qq{\nUsage: $0 -mode initialise -file <filename> -dbname <db> }
          .qq{-dbhost <host> -dbuser <user> -dbport <port> -dbpass <pass>\n};
        exit(0);
    }
    $threads = 17;
    $initial_populate = 1;
    $sbmid = 0;
    $user = 'admin';
    $uid = 0;
} else {
    if (!$sbm_file || !$db    || !$dbhost || !$dbuser|| !$dbpass || !$dbport ||
        !$sbmid    || !$user  || !$uid) {
        print qq{\nUsage: $0 -file <filename> -dbname <db> 
    -dbhost <host> -dbuser <user> -dbport <port> -dbpass <pass>
    -sbmid <sbmid> -capuser <capuser> -capuid <uid>\n} ;
    exit(0);
    }
    $threads||=1;

    my @time = localtime;
    my @mon = qw/Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec/;
    $time[5] -= 100;
    my $tstamp = $time[3].$mon[$time[4]].$time[5];
    my $logfile = 'gff_'.$tstamp.'.log';
    my $appender2 = Log::Log4perl::Appender->new("Log::Dispatch::File",filename=>$logfile,mode=>"append");
    ##my $layout = Log::Log4perl::Layout::PatternLayout->new("%d %p> %F{1}:%L %M - #%m%n");
    $log4->add_appender($appender2);
}

my $dbstr = qq{DBI:mysql:database=$db;host=$dbhost;port=$dbport;}
              .qq{user=$dbuser;password=$dbpass};
my $dbhashref = {PrintError => 1, RaiseError => 1};
my $dbi = DBI->connect($dbstr, $dbhashref) or die;

$log4->info(qq{processing $sbm_file from $user/$uid for sbm $sbmid});

#### removing multiple file-handling - all will be done by wrapper?!? - can't be arsed to chage data structure?!?
#@file_list_arrayofhashref = (+{file => $sbm_file, id => 0, user => 'admin', uid => 0}); 
my $sbm = +{file => $sbm_file, id => $sbmid, user => $user, uid => $uid}; # + is totally unecessary...

if ($initial_populate) {
    print '=' x ($interval);
    #print qq{\n};
    $log4->info('truncating feat table');
    my $clean_sth = $dbi->prepare('truncate table feat');
    $clean_sth->execute() or die;
    $log4->info('truncating feat_rel table');
    $clean_sth = $dbi->prepare('truncate table feat_rel');
    $clean_sth->execute() or die;
    $log4->info('truncating canon_crc64 table');
    $clean_sth = $dbi->prepare('truncate table canon_crc64');
    $clean_sth->execute() or die;
} else {

    # put in a check for whether a submission id has been processed before?!? i.e.
    # select where sbmid =1 and status = old...

    my $sth = $dbi->prepare('select * from feat where sbm_id = '.$sbmid
      .' and status = "old"');
    $sth->execute() or die;
    #/ no real need for fetchall - any result is a response?!?
    if (@{$sth->fetchall_arrayref()}) {

        print qq{\n[ERROR] There are already entries within the data base with submission ids }
          .qq{corresponding to this one?!?\n};
        exit(1);
        #exit(0);
    }
}    

# for my $sbm (@file_list_arrayofhashref) { #while (my $sbm = $sth->fetchrow_hashref) {
#my $sbm_file = $upload_dir.$sbm->{file};

#y/ possibly move to sys::mmap?!? seems to be issue with opening file so intead slurp the fucker?!?
#y/ and fragment - the files are never gonna be all that big?!? that ain't the hold up!?!
#={{{ sys::mmap stuff
    ## use Sys::Mmap;

    ##new Mmap $str, 8192, 'structtest2.pl' or die $!;
    ##new Mmap $var, 8192 or die $!;

    ##mmap($foo, 0, PROT_READ, MAP_SHARED, FILEHANDLE) or die "mmap: $!";
    ##@tags = $foo =~ /<(.*?)>/g;
    ##munmap($foo) or die "munmap: $!";
    
    ##mmap($bar, 8192, PROT_READ|PROT_WRITE, MAP_SHARED, FILEHANDLE);
    ##substr($bar, 1024, 11) = "Hello world";

    ##mmap($baz, 8192, PROT_READ|PROT_WRITE, MAP_SHARED|MAP_ANON, STDOUT);

    ##$addr = mmap($baz, 8192, PROT_READ|PROT_WRITE, MAP_SHARED|MAP_ANON, STDOUT);
    ##Sys::Mmap::hardwire($qux, $addr, 8192);
#=}}}

my $filebytes = -s $sbm_file;
#my $hm = `ls -al $sbm_file`;
my $fh;

# allow for empty files - i.e. 0
if (!defined $filebytes) {
    $log4->error(qq{couldn't check size of file $sbm_file});
    $log4->info(qq{exiting});
    exit(1);
} elsif (!$filebytes) {
    $log4->error(qq{file $sbm_file is empty});
    $log4->info(qq{exiting});
    exit(1);
} elsif (!open($fh, '<', $sbm_file)) {
    $log4->error(qq{couldn't open file $sbm_file for reading});
    $log4->info(qq{exiting});
    exit(1);
}

#$filebytes || die qq{\n$capproc couldn't check size of file $sbm_file\n$capproc exiting};

my $pieces = $threads > 1 ? $threads - 1 : 0;
my $chunk = $filebytes/($threads);
#my $chunk = length($map)/($threads);

#while ($map =~ /(\n)/g) {  #print "Word is '$1', ends at position ", pos $map, "\n";  #}





my @offsets;
for (0..$pieces-1) { #while (tell($fh) < length($map)) {
    seek $fh,$chunk,1; # we set new position to where we are + chunk size and not to absolute position with ,,0
    my $g = $fh->getline;
    push @offsets, tell($fh); # print q{\nhence we know we are at the start of a line at }.tell($fh);
}

close $fh;
push @offsets, $filebytes;
unshift @offsets, 0;

#print qq{\noffset: $_} for (@offsets);

map_file(my $map, $sbm_file);
my @ths;
$log4->info(qq{creating $threads threads});
for (0..$#offsets-1) {
    my $offset = $offsets[$_];
    my $length = $offsets[$_+1]-$offset;
    my $bit = substr $map, $offset, $length;
    push @ths, threads->create('_start_processing_thread', $bit, $sbm, $_);
}

for (@ths) {
    $_->join();
}

$log4->info(q{joined threads});

#b now handled by wrapper script?!?
#if ($c==0) {
#    print qq{\nINFO: No new files to import - exiting!\n};
#    exit(0);
#}

@permitted_biotypes = map { qq{'$_'} } @permitted_biotypes;
my $permitted_biotypes = join(q{,}, @permitted_biotypes);

$log4->info(q{processing crc64 values for gene models});
my $sth = $dbi->prepare(qq{select 
  f.sbm_id, f.name, f.f_id, f.seqid, f.strand, f.start, f.end 
  from feat f where
  f.status = 'new' and
  f.type = 'gene'});
$sth->execute() or die;

my $keep_feature_list = [];
my $keep_feature_rel_list = [];

#y we don't give a crap about the details here we just want them in a list to thread the thing?!?
#my @genes_to_check = map { $_->{f_id} } 

#while (my $gene = $sth->fetchrow_hashref) {

my $out = $sth->fetchall_hashref('f_id');
my @list = (values %{$out});

my $inty = int((scalar @list)/$threads);
my $pos = 0;
@ths = ();

$log4->info(qq{creating $threads threads});
for my $t (0..$threads-1) {
    my $start = $t == 0 ? 0 : $pos+1;
    my $end = $t == $threads-1 ? $#list :$pos + $inty;
    $pos = $end;
    my @sub = @list[$start..$end];
    push @ths, threads->create('_thread_crc', $t, @sub);
    #$pos += $end;
}

for (@ths) {
    $_->join();
}

$log4->info(q{joined threads});

$log4->info(qq{adding crc64 index for canonical gene set});

exit if ($insert_only);

if ($initial_populate) {
    $log4->info(qq{ re-truncating feat table});
    my $clean_sth = $dbi->prepare('truncate table feat');
    $clean_sth->execute() or die;
    $log4->info(qq{ re-truncating feat_rel table});
    $clean_sth = $dbi->prepare('truncate table feat_rel');
    $clean_sth->execute() or die;
    $log4->info(qq{ database is ready for cap submissions});
    my $sth = $dbi->prepare('insert into status (status) values ("ready"');
    $sth->execute() || die;
} else {

    #r rather than having status on feat_rel too and/or keeping track of f_id to bin do table join for wiping!?!

    # we ban top-level features that aren't genes. this means that we end up with
    # orphan offspring features - to remove these clearly a table join won't work
    # thus we either have to relax the gene condition, wipe orphan features - way
    # inefficient cos of all the table joins or put a status column in feat_rel too?!?

    # serious table joins:  select distinct f1.name, f1.f_id,f1.type from feat f1 left join (feat f2, feat_rel fr) on (f1.name=f2.name and f2.name=fr.parent_name) where f1.name is not null;
    # select distinct f1.name, f1.f_id,f1.type from feat f1 left join (feat f2, feat_rel fr) on (f1.name=f2.name and f2.name=fr.parent_name) where f1.name is null;

    $log4->info(qq{protected $modified new/modified gene models (file contained $features permited feature lines)});

    $log4->info('removing all unprotected features/relationships');
    my $delete_others_sth = $dbi->prepare(qq{delete fr.* 
      from feat_rel fr 
      where fr.status = 'new'}
    );
    $delete_others_sth->execute() || exit(1);
    $delete_others_sth = $dbi->prepare(qq{delete f.*
      from feat f
      where f.status = 'new'}
    );
    $delete_others_sth->execute() || exit(1);

    #y update db
    my $updatedb = $dbi->prepare("update files set 
        feature_no = $features, 
        new_mod_genes = $modified 
        where sbm_id = $sbmid"
    );
    $updatedb->execute() || exit(1);
}

$sth->finish();

#y if we got here we tell wrapper all okay?!?
exit(0);

sub print_number {
    my ($str) = @_;
    system qq{echo -n $str};
    return;
}

sub print_clean {
    my $l = shift;
    print qq{\x08} x $l;
    return;
}

#mysql> alter table files modify column status enum('new','old','prob');
sub schema {

print <<BLAH

DROP TABLE IF EXISTS `status`;
CREATE TABLE `status` (
  `status_id` int(11) NOT NULL AUTO_INCREMENT,
  `ts` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `status` enum ('ready', 'looping', 'error') NOT NULL,
  PRIMARY KEY (`status_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `meta`;
CREATE TABLE `meta` (
  `meta_id` int(11) NOT NULL AUTO_INCREMENT,
  `meta_key` varchar(15) NOT NULL,
  `meta_value` varchar(15) NOT NULL,
  PRIMARY KEY (`meta_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
insert into meta (meta_key,meta_value) values ('status','ready');

DROP TABLE IF EXISTS `files`;
CREATE TABLE `files` (
  `sbm_id` int(11) NOT NULL AUTO_INCREMENT,
  `user` varchar(15) NOT NULL,
  `uid` int(11) NOT NULL,
  `mail` varchar(40) NOT NULL,
  `description` varchar(50) DEFAULT NULL,
  `md5` varchar(40) NOT NULL,
  `size` int(11) NOT NULL,
  `hostname` varchar(40) NOT NULL,
  `species` enum('Phumanus','Gmorsitans'),
  `file` varchar(50) NOT NULL,
  `type` enum('gff3','fasta','xls') not null,
  `ts` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `status` enum('new','old', 'prob', 'not_gff', 'empty') NOT NULL,
  `feature_no` int(11) DEFAULT NULL,
  `new_mod_genes` int(11) DEFAULT NULL,
  `tmp_gff` varchar(40) default NULL,
  CONSTRAINT uniq_file UNIQUE (md5),
  PRIMARY KEY (`sbm_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

#### amazingly some genes were breaking the system as their names were >40 chars?!?
DROP TABLE IF EXISTS `feat`;
# no idea of exact score encoding?!?
CREATE TABLE `feat` (
  `f_id` int(11) NOT NULL AUTO_INCREMENT,
  # make an enum?!?
  `type` varchar(10) DEFAULT NULL,
  `seqid` varchar($name_length) NOT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  `strand` enum('.','+','-') NOT NULL,
  `source` varchar(10) NOT NULL,
  `score` varchar(10) NOT NULL,
  `phase` enum('.','0','1','2') NOT NULL,
  `annot` varchar(256) NOT NULL,
  `name` varchar($name_length) DEFAULT NULL,
  `sbm_id` int(11) NOT NULL,
  `status` enum('new','old') NOT NULL,
  CONSTRAINT uniq_in_sub UNIQUE (name,sbm_id),
  PRIMARY KEY (`f_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
#) ENGINE=InnoDB AUTO_INCREMENT=501 DEFAULT CHARSET=latin1;

# UNIQUE KEY `crc64` (`crc64`)
# `canonical` int(11) DEFAULT NULL,
# `crc64` varchar(60) DEFAULT NULL,

DROP TABLE IF EXISTS `canon_crc64`;
CREATE TABLE `canon_crc64` (
  `crc64` varchar(50) NOT NULL,
  PRIMARY KEY (`crc64`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
    
DROP TABLE IF EXISTS `feat_rel`;
CREATE TABLE `feat_rel` (
  `r_id` int(11) NOT NULL AUTO_INCREMENT,
  `f_id` int(11) NOT NULL,
  `parent_name` varchar($name_length) NOT NULL,
  `status` enum('new','old') NOT NULL,
  PRIMARY KEY (`r_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

BLAH

}

sub _start_processing_thread {

    my ($bits, $sbm, $id) = @_;

    #y clearly each thread needs its own db connection so we scope the name...
    my $dbi = DBI->connect($dbstr, $dbhashref);

    my $sbm_file = $bits;
    my $sbm_id = $sbm->{id};
    my $sbm_user = $sbm->{user};
    my $sbm_uid = $sbm->{uid};
    my $linen = 0;
    my @lines = $bits =~ /([^\n]+)\n/g;

    my $l = scalar @lines;
    my $int_lines = int(($l)/$interval) != 0 ? int(($l)/$interval) : 1;

    my $local_feat = 0;
    LINE:
    for my $line (@lines) {
    #while (my $line = <$fh>) {
    #while ($bits =~ /([^\n]+\n)/g) { 

        #y just print in one thread - can't be bothered with more shared vars...
        # print_number('=') if ($int_lines != 0 && $linen % $int_lines == 0 && $id == 0);
        $linen++;
        print_number(sprintf("%02s%%",int(100*$linen/$l))) 
        if ($int_lines != 0 && $linen % $int_lines == 0 && $id == 0);
        &print_clean(3);
            
        next LINE if ($line =~ /^[\x20\x09]*$/x); # ignore empty lines
        next LINE if ($line =~ /^[\x20\x09]*\x23/x); # ignore comment lines

        if(my @vals = $line =~/^
          ([^\t]+)\t   # seqid
          ([^\t]+)\t          # source
          ([^\t]+)\t   # type
          ([^\t]+)\t   # start
          ([^\t]+)\t   # end
          ([^\t]+)\t          # score
          ([^\t]+)\t   # strand
          ([^\t])+\t          # phase
          ([^\t]+)     # attribs
          $/xg) { #/ should compile it with /o

            my $type = $vals[2];

            #/ parentless utrs etc.?!?
            next LINE if ($type =~ /^($prohibited_types)$/);

            $local_feat++;
            #////// not really necessary, just stops wasting time?!?
            # next LINE if ($vals[1] !~ /^(gene|CDS|exon|$permitted_biotypes)/);
            # next LINE if ($vals[1] eq 'protein' || $vals[1] eq 'polypeptide');

            my $comment = $vals[8]; #y each new match resets $2/$4...

            #y handle either comma-delimited parent names and multiple parent declarations
            my @parents = ($comment =~ /Parent=\s?(\S+?)\s?;/g);
            if ((scalar @parents == 1) && ($parents[0] =~ /,/)) {
                @parents = split(q{,}, $parents[0]);              
            }

            my $annotstr;
            if ($line =~ /(.+?)([^\t]+)$/) {
                $line = $1;
                $annotstr = $2;
            }

            my $id_str = q{};
            if ($comment =~ /ID=(\S+?)\s?\;/) {
                $id_str = $1;
                die qq{\nNames longer than 256 char are not supported atm?} if (length($id_str) > 256);
                $id_str = qq{, name = '$1' };
            }

            # absolutely refuse any top level features that aren't 'gene' can't allow 
            # this to slip as it will make gff from multiple source too messy to deal with

            if (
                #y ban all but gene and sometimes pseudogene from having id without parents - i.e. top-level
                #y or floating features?!?
                ($id_str && (scalar @parents == 0))
                && ($type ne 'gene' && $type ne 'pseudogene') 
                # (($type eq 'pseudogene') && !$wtsi))
            ) {
                next LINE if ($type eq 'protein' || $type eq 'polypeptide');
                #/ just cos artemis seems to use pseudogene/pseudogenic exon etc.,...
                # next LINE if ($wtsi && $type eq 'pseudogene');
                print qq{\nonly 'gene' features are permitted as top-level features!};
                print qq{\nid: $id_str\nparents: @parents\ntype:$type\n};
                die;
            }

            my $store_feat = qq{INSERT INTO feat SET seqid = ?, start = ?, end = ?,  }  
            .qq{strand = ?, phase = ?, annot = ?, type = ?, sbm_id = ?, status = 'new', }
            .qq{source = ?, score = ? }
            #/ can't be bothered with conditional bind_param
            .$id_str;

            $sth = $dbi->prepare($store_feat);
            $sth->bind_param( 1, $vals[0],  SQL_VARCHAR );
            $sth->bind_param( 2, $vals[3],  SQL_INTEGER );
            $sth->bind_param( 3, $vals[4],  SQL_INTEGER );
            $sth->bind_param( 4, $vals[6],  SQL_CHAR );
            $sth->bind_param( 5, $vals[5],     SQL_VARCHAR );
            #$sth->bind_param( 5, $line,     SQL_VARCHAR );
            $sth->bind_param( 6, $annotstr, SQL_VARCHAR );
            $sth->bind_param( 7, $type,  SQL_VARCHAR );
            $sth->bind_param( 8, $sbm_id,   SQL_INTEGER );
            $sth->bind_param( 9, $vals[1],  SQL_VARCHAR );
            $sth->bind_param( 10, $vals[7],  SQL_VARCHAR );
            $sth->execute();

            # can just use last_insert_id() in insert statemet - but we will reuse
            my $fid = $sth->{'mysql_insertid'}; 
            
            for my $p (@parents) {
                #### must have sbm-id correlate and the combo of id and sbm-id must be unique
                my $store_feat_rel = qq{INSERT INTO feat_rel SET f_id = ?, parent_name = ?, status = 'new' };
                $sth = $dbi->prepare($store_feat_rel) or die;
                $sth->bind_param( 1, $fid, SQL_VARCHAR );
                $sth->bind_param( 2, $p,   SQL_VARCHAR );
                $sth->execute() or die;
            }

        } else { 
            die qq{\n* non-legal gff line at $. in file $sbm_file\n}; 
        }
    }

    $features += $local_feat;

    return;
}

sub _thread_crc {

    my $id = shift;
    my @sub_genes = @_;
    my $dbi = DBI->connect($dbstr, $dbhashref);

    #={{{ some more statement prep
    # ONLY transcripts should have parents with gene names e.g. utrs etc., have
    # transcript parents so can realx this constraint - otherwise need to know EVERY type!?!
    my $gene_status1_update_sth = $dbi->prepare(qq{update feat set status = 'old' where f_id = ?});
    my $gene_status2_update_sth = $dbi->prepare(qq{update feat_rel set status = 'old' where r_id = ?});
    my $crc_check_sth = $dbi->prepare(qq{select * from canon_crc64 where crc64 = ?});
    my $crc_insert_sth = $dbi->prepare(qq{insert into canon_crc64 set crc64 = ? });
    my $mrna_sth = $dbi->prepare(qq{select 
      f.name, f.seqid, f.strand, f.start, f.end, f.f_id, fr.r_id
      from feat f, feat_rel fr
      where fr.f_id=f.f_id 
      and sbm_id = ? 
      and fr.parent_name = ? 
    });
    # and type in ($permitted_biotypes)});
    # and type = 'mRNA'});
    my $termini_sth = $dbi->prepare(qq{select 
      seqid, f.strand, f.start,f.end,f.type,f.f_id, fr.r_id
      from feat f,feat_rel fr
      where f.f_id=fr.f_id and 
      f.sbm_id = ? and 
      f.type = ? and 
      fr.parent_name = ?
    });
    #=}}}

    my $c = 0;
    my $l = scalar @sub_genes;
    my $int_lines = int(($l)/$interval) != 0 ? int(($l)/$interval) : 1;

    for my $gene (@sub_genes) {
    # while (my $gene = $sth->fetchrow_hashref) {

        $c++;
        print_number(sprintf("%02s%%",int(100*$c/$l))) 
          if ($int_lines != 0 && $c % $int_lines == 0 && $id == 0);
        &print_clean(3);

        my @gene_feat_list;
        my @gene_feat_rel_list;
        my $g_fid = $gene->{f_id};
        my $g_name = $gene->{name};
        my $g_sbmid = $gene->{sbm_id};
        my $g_string = $cd1.$gene->{seqid}.$cd.$gene->{strand}.$cd.$gene->{start}.$cd.$gene->{end};

        push @gene_feat_list, $g_fid;

        #y same for exons and cds so do it now...
        $termini_sth->bind_param( 1, $g_sbmid, SQL_INTEGER );

        $mrna_sth->bind_param( 1, $g_sbmid, SQL_INTEGER );
        $mrna_sth->bind_param( 2, $g_name,  SQL_VARCHAR );
        $mrna_sth->execute() or die;

        my @mrna_strings;
        while (my $mrna = $mrna_sth->fetchrow_arrayref) {
        
            my $t_name = $mrna->[0]; # sbm_id will be same as gene's
            push @gene_feat_list, $mrna->[5];
            push @gene_feat_rel_list, $mrna->[6];

            my $type = 'exon';
            $termini_sth->bind_param(     2,      $type,       SQL_VARCHAR );
            $termini_sth->bind_param(     3,      $t_name,       SQL_VARCHAR );
            $termini_sth->execute();
            my $exons = $termini_sth->fetchall_arrayref() or die;

            $type = 'CDS';
            $termini_sth->bind_param(     2,      $type,       SQL_VARCHAR );
            $termini_sth->execute();
            my $cdss = $termini_sth->fetchall_arrayref() or die;

            my @exon_strings;
            my @cds_strings;
            for my $p (@{$exons}) {

                ##### why the fuck is type going in?!?
                push @exon_strings, $cd3.$p->[0].$cd.$p->[1].$cd.$p->[2].$cd.$p->[3];
                #push @exon_strings, $p->[0].$cd.$p->[1].$cd.$p->[2].$cd.$p->[3];
                #push @exon_strings, $p->[0].$cd.$p->[1].$cd.$p->[2].$cd.$p->[3].$cd.$p->[4];
                push @gene_feat_list, $p->[5];
                push @gene_feat_rel_list, $p->[6];
            }
            for my $p (@{$cdss}) {
                push @cds_strings, $cd4.$p->[0].$cd.$p->[1].$cd.$p->[2].$cd.$p->[3];
                #push @cds_strings, $p->[0].$cd.$p->[1].$cd.$p->[2].$cd.$p->[3];
                #push @cds_strings, $p->[0].$cd.$p->[1].$cd.$p->[2].$cd.$p->[3].$cd.$p->[4];
                push @gene_feat_list, $p->[5];
                push @gene_feat_rel_list, $p->[6];
            }
        
            @exon_strings = sort {$a cmp $b} @exon_strings;
            @cds_strings = sort {$a cmp $b} @cds_strings;
            my $es = join(q{},@exon_strings);
            my $cs = join(q{},@cds_strings);

            push @mrna_strings, $cd2.$mrna->[1].$cd.$mrna->[2].$cd.$mrna->[3].$es.$cs;
        }

        @mrna_strings = sort {$a cmp $b} @mrna_strings;
        my $ms = join(q{},@mrna_strings);
        $g_string .= $ms;

        $termini_sth->finish();
        $mrna_sth->finish();

        my $g_crc64 = `$crc64 $g_string`;

        if ($crc_debug && $g_name =~ $crc_debug) {
            $g_string =~ s/1@/\n[gene]\t/g;
            $g_string =~ s/2@/\n\t[tran]\t/g;
            $g_string =~ s/3@/\n\t\t[exon]\t/g;
            $g_string =~ s/4@/\n\t\t[cdss]\t/g;
            print qq{\n[start] search-name: $crc_debug stored-name: $g_name${g_string}\n[end] -> $g_crc64\n};
        }

        if($initial_populate) {
            $crc_insert_sth->execute($g_crc64) or die;
        } else {
            $crc_check_sth->execute($g_crc64) or die;

            #y we give safe status to features and wipe all others to stop everything getting messy?!?
            if (@{$crc_check_sth->fetchall_arrayref()}) { # unique constraint on table so can only have one return per sbm!?!
                @gene_feat_list = ();
            } else {

                # If a container object, such as a hash or array, is locked, all the
                # elements of that container are not locked. For example, if a thread does
                # a lock @a, any other thread doing a lock($a[12]) won't block.

                for my $fid (@gene_feat_list) {
                    $gene_status1_update_sth->execute($fid) or die;
                }    
                for my $rid (@gene_feat_rel_list) {
                    $gene_status2_update_sth->execute($rid) or die;
                }

                lock($modified);
                $modified++;
            }
        }

    }
    return;
}    
