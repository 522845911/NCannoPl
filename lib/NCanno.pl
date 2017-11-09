#!/usr/bin/env perl
#$ -S /usr/bin/perl
#$ -cwd
# Copyright (c) GenePlus 2015
# Program       : NCanno.pl
# Author        : Liu Tao
# Program Date  : 2015-07-08
# Modifier      : 
# Last Modified : 
# Description   : annotate vcf by BedAnno, and database API, give a well formatted serious mutants, 
#                 grouped by transcript id, together with the protein id
# Dependency    : This script depend on the bgzip,tabix and Faidx module
#                 to be available, and the whole genome fasta to be read from 
#                 db/anno/aln_db/hg19/ named "hg19_chM.fa", which change 
#                 the original chrM of hg19 to chrM_NC_012920.1.
#                 Also the BedAnno annotation databases are needed to annotate.

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Spec;
use File::Basename;
use Cwd 'abs_path';
use Carp;
use threads;
use threads::shared;
use Thread::Queue;
use Time::HiRes qw(gettimeofday tv_interval);
use IO::Uncompress::Gunzip qw($GunzipError);

our $NOAHCARE;
BEGIN {
    $NOAHCARE = dirname(dirname(dirname(abs_path($0)))) . '/program/NoahCare/';
}

my $VERSION = '1.11';

use lib $NOAHCARE.'/lib';
use BedAnno;

my $beda;
# localize SIGNAL handler
$SIG{'INT'}     = \&quit;
$SIG{'QUIT'}    = \&quit;
$SIG{'FPE'}     = \&quit;
$SIG{'ABRT'}    = \&quit;
$SIG{'TERM'}    = \&quit;


# database file default name configuration
my $hg19_fa   = $NOAHCARE.'/db/anno/aln_db/hg19/hg19_chM.fa.rz';
my $db        = $NOAHCARE.'/db/anno/annodb/ncbi_anno_rel104_db.bed.gz';
my $trDB      = $NOAHCARE.'/db/anno/annodb/ncbi_anno_rel104_trSeq.fas';
my $cytoband  = $NOAHCARE.'/db/anno/cytoBand/cytoBand_hg19_grch37.txt.gz';
my $pfam      = $NOAHCARE.'/db/anno/pfam/Pfam-A-ncbi_2012-12-21.bed.gz';
my $preds     = $NOAHCARE.'/db/anno/predictions/predictDB_for_anno104.tab.gz';
my $condel    = $NOAHCARE.'/config/Condel/';
my $phyloP    = $NOAHCARE.'/db/anno/phyloP/phyloP_all3class_combin_2013-09-25.bed.gz';
my $dbsnp     = $NOAHCARE.'/db/anno/dbsnp/snp147.bed.gz';
my $esp6500   = $NOAHCARE.'/db/anno/NHLBI/ESP6500SI-V2-SSA137.NHLBI.bed.rmanchor.uniq.gz';
my $tgp       = $NOAHCARE.'/db/anno/tgp/tgpPhase3_vcfDB/tgp_phase3_small_vars.vcf.gz';
my $repeat_db = $NOAHCARE.'/db/anno/rep_dup_db/dup_trf_rmsk.bed.gz';
my $omim_data = $NOAHCARE.'/db/anno/omim/anno_combin.tsv.gz';
my $cgd_data  = $NOAHCARE.'/db/anno/CGD/CGD.txt.gz';
my $cosmic    = $NOAHCARE.'/db/anno/cosmic/COSMICv80_GRCh37_20170213.bed.gz';
my $exac_data = $NOAHCARE.'/db/anno/exac/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz';
my $gnomAD_db = $NOAHCARE.'/db/anno/gnomAD/gnomad.exomes.r2.0.1.sites.vcf.gz';
my $hgmd      = $NOAHCARE.'/db/anno/hgmd/localhgmd.db.gz';
my $verList   = $NOAHCARE.'/config/db_version.list';
my $msqc_info = $NOAHCARE.'/db/analyze/snv/MSQC_info.list';
my $clinVar   = $NOAHCARE.'/db/anno/clinVarDB/clinVarDB.vcf.gz';
my $vusVar    = $NOAHCARE.'/db/anno/clinVarDB/vusVarDB.vcf.gz';

my %DB_VERSION = (
    REFERENCE_BUILD => 'GRCh37',
    ANNOTATION_BUILD => 'NCBI Annotation Release 104',
    ENSDB_BUILD => '73',
    DBSNP_BUILD => '147',
    COSMIC_BUILD => '80',
    TGP_BUILD => 'phase3_v1',
    CGD_Date => '2013-12-30',
    OMIM_Date => '2015-01-19',
    PFAM_Date => '2012-12-21',
    LOCALHGMD_YEAR => '2012',
    ESP6500_VERSION => '2',
    ExAC_VERSION => '0.3.1',
    GAD_VERSION => '2.0.1',
    UCSC_HG19_Date => '2013-12-06',
    CLIVAR_Date => '2016-11-28',
);

my $finished; # don't define it any where until finished.

# check necessary dependencies
check_dep( $db, $trDB, $hg19_fa );

my (
    $help,         $outfile,       $sample_list, $quiet,
    $headrule,     $type,          $offline,     $number_threads,
    $buffer_cache, $write_version, $debug,       $timer,
    $in_msqc,      $msqc_rst
);

# running mode
$type        = "vcf";
my $NUM_THR     = 2;
my $BUFFER_SIZE = 5000;
my $MAX_REF_LEN = 200;

# filter rule thresholds
my $AD_UP_THRESHOLD = 8;  # used to decide PASS
my $AD_DN_THRESHOLD = 2;  # used to decide FAIL
my $PL_THRESHOLD    = 3;  # used to decide PLTAG (FAIL)
my $PL_UP_THRESHOLD = 10; # used to decied PLTAG (DUBIOUS)

my %deleterious_func;
@deleterious_func{qw(
    cds-loss init-loss stop-gain stop-loss
    frameshift knockout nonsense
)} = (1) x 7;

my %possible_deleterious_func;
@possible_deleterious_func{qw(
    missense stop-retained splice splice-3 splice-5 span abnormal-intron
)} = (1) x 7;

my %autoInterp_freq_threshold = (
    TGP_AF    => 0.01,
    ESP6500AF => 0.01,
    PVFD_AF   => 0.01,
    ExAC_AF   => 0.01,
    GAD_AF    => 0.01,
    PanelAF   => 0.05,
);

my %kickoutExcel_function;
@kickoutExcel_function{qw(
    utr-3 utr-5 promoter intron annotation-fail unknown-no-call unknown no-change .
)}= (1) x 9;

my %AutoInterpWords = (
    Certain => 1,
    "Likely Deleterious" => 2,
    "VUS" => 3,
    "Likely Benign" => 4,
    Benign => 5,
    Unknown => 6,
);

my @xpar_reg = (
    [ 60000, 2699520 ],
    [ 154931043, 155260560 ]
);

my @ypar_reg = (
    [ 10000, 2649520 ],
    [ 59034049, 59363566 ]
);

my %VarNameList;

GetOptions(
    "type|t=s"       => \$type,
    "numthreads|n=i" => \$number_threads,
    "buffersize|b=i" => \$buffer_cache,
    "samplist|s=s"   => \$sample_list,
    "outfile|o=s"    => \$outfile,
    "headrule|r=s"   => \$headrule,
    "inmsqc|i=s"     => \$in_msqc,
    "qcresult|c=s"   => \$msqc_rst,
    "offline|f"      => \$offline,
    "quiet|q"        => \$quiet,
    "debug|d"        => \$debug,
    "timer|m"        => \$timer,
    "writever|w"     => \$write_version,
    "help|h"         => \$help
);

require Vcf if ($type ne "tsv");
pod2usage( -verbose => 2, -exitval => 0 ) if ($help);
(
         !@ARGV
      or !-r $ARGV[0]
      or $ARGV[0] =~ /\.vcf(\.gz)?$/i
      or $ARGV[0] =~ /\.tsv(\.bz2)?$/i
) and pod2usage( -verbose => 1, -exitval => 1 );

our %opts; # BedAnno options
# some db is not available online currently, 
# these local version database APIs have been 
# integrated into BedAnno
check_dep(
    $cytoband,  $pfam,   $preds,  $condel, $repeat_db, $esp6500,
    $omim_data, $phyloP, $cosmic, $dbsnp,  $exac_data, $gnomAD_db,
    $tgp
);
$opts{cytoBand}   = $cytoband;
$opts{pfam}       = $pfam;
$opts{prediction} = $preds;
$opts{condel}     = $condel;
$opts{esp6500}    = $esp6500;
$opts{phyloP}     = $phyloP;
$opts{cosmic}     = $cosmic;
$opts{rmsk}       = $repeat_db;
$opts{dbSNP}      = $dbsnp;
$opts{genome}     = $hg19_fa;
$opts{exac}       = $exac_data;
$opts{gnomAD}     = $gnomAD_db;
$opts{tgp}        = $tgp;

# Omim and CGD is outof the scope of BedAnno
check_dep($cgd_data, $hgmd, $clinVar, $vusVar);

my $OutFp = \*STDOUT;
if ($outfile) {
    open OUT, ">", $outfile or confess "$outfile: $!\n";
    $OutFp = \*OUT;
}

if (defined $msqc_rst) {
    open QCSTAT, ">", $msqc_rst or confess "$msqc_rst: $!\n";
}

if ($headrule) {
    open HDR, ">", $headrule or confess "$headrule: $!\n";
}

# read the config items
our %config;
require "$ARGV[0]";

# Annotation rules
# available pop: EAS SAS EUR AMR AFR
$config{TGP_POP}       = "EAS"   if ( !exists $config{TGP_POP} );

# available pop: EAS,SAS,AMR,AFR,EUR,FIN,NFE,OTH
$config{ExAC_POP}      = "EAS"   if ( !exists $config{ExAC_POP} );
$config{Panel_Control} = "."     if ( !exists $config{Panel_Control} );
$config{Panel_ID}      = "."     if ( !exists $config{Panel_ID} );
$config{FLANK_LEN}     = 2       if ( !exists $config{FLANK_LEN} );
$config{intron_edge_remain} = 0  if ( !exists $config{intron_edge_remain} );

$verList = $config{DBVER_LIST} if ( exists $config{DBVER_LIST} );
$msqc_info = $config{MSQC_INFO} if ( exists $config{MSQC_INFO} );
$NUM_THR = $config{ANNO_THR_NUM} if ( exists $config{ANNO_THR_NUM} );
$BUFFER_SIZE = $config{ANNO_BUFFER_SIZE}
  if ( exists $config{ANNO_BUFFER_SIZE} );

# command line option have the highest priority
$NUM_THR = $number_threads if (defined $number_threads);
$BUFFER_SIZE = $buffer_cache if (defined $buffer_cache);

$NUM_THR = 2 if ($NUM_THR < 2);
$BUFFER_SIZE = 100 if ($BUFFER_SIZE < 100);
my $SINGLE_TIMEOUT = 60 * ($NUM_THR - 1);

if ( exists $config{autoInterp_freq_threshold} ) {
    %autoInterp_freq_threshold = %{$config{autoInterp_freq_threshold}};
}

if ( exists $config{kickoutExcel_function} ) {
    %kickoutExcel_function = %{$config{kickoutExcel_function}};
}

if ( exists $config{possible_deleterious} ) {
    @possible_deleterious_func{ @{ $config{possible_deleterious} } } =
      (1) x ( scalar @{ $config{possible_deleterious} } );
}

if ( exists $config{VarNameList} ) {
    open (VNL, $config{VarNameList}) or die "Error: [$config{VarNameList}] $!\n";
    while (<VNL>) {
        next if (/^#/ or /^\s*$/);
        s/\s+$//;
        my ($g, $v, $n) = split(/\t/,$_,4);
        $VarNameList{$g}{$v} = $n;
    }
    close VNL;
}

if ( exists $config{CoreFuncRegion} ) {
    check_dep( $config{CoreFuncRegion} );
}


our $GeneTestCode;
if (exists $config{GeneTestCode}) {
    $GeneTestCode = read_TestCode($config{GeneTestCode});
}

# init Homo reference var list
my ($HomoRef_Var, $HRVanno_opt);
if (exists $config{HomoRef_Var}) {
    $HomoRef_Var = read_HomoRef_var($config{HomoRef_Var});
}

# read and init MS QC info
my ($QCsites, $msQCrst, $ngsQCrst, $nQCs);
if (defined $in_msqc) {
    $QCsites = read_MSQC_info($msqc_info);
    $msQCrst = read_inmsqc($in_msqc);
    my @msqc_samples = keys %$msQCrst;

    my @sids = ();
    foreach my $chr (keys %$QCsites) {
        push (@sids, sort keys %{$QCsites->{$chr}});
    }
    my %qcsids = map { $_ => 1 } @sids;
    $nQCs = scalar keys %qcsids;

    unless ( 0 < @msqc_samples
        and $nQCs == scalar( @{ $$msQCrst{ $msqc_samples[0] } } ) )
    {
        die join( "\n",
            "Error: inconsistent - QC sites info file and the MS QC result",
            "       QC sites info file : $msqc_info",
            "       MS QC result file  : $in_msqc\n" );
    }
}

our $sample_info = {};
if ($sample_list) {
    $sample_info = read_sample_list($sample_list);
}

if (exists $config{"genes"}) {
    check_dep( $config{"genes"} );
    $opts{genes} = $config{"genes"};
}
if (exists $config{"trans"}) {
    check_dep( $config{"trans"} );
    $opts{trans} = $config{"trans"};
}

print STDERR "== Annotation Start at " . localtime() . " ==\n";

# shift the first args of config file
my $CONFIG_FILE = shift(@ARGV);

my %var_count = ();
my %no_call_var_count = ();
my %sex_varcount = ();
my %ti = (); # A <=> G, C <=> T
my %tv = (); # others
my %annotation_count = ();
my %pre_total_AD = ();

my @out_header_keys;
# print the header line, assign out_header_keys
# and output the format of variation sheet if needed
print_header();

our @anno_cache : shared;
@anno_cache = ();

# The input variation file must be carefully sorted
our $cached_count : shared;
$cached_count  = 0;

# Annotation Queue
my $anno_q = Thread::Queue->new();
my @anno_thrs = map {threads->new( \&anno_threads )} (1 .. ($NUM_THR - 1));

my %prev_var;
my $enqueue_count = 0;
my $read_record_time = 0;
my $main_record_count = 0;
my $total_record_count = 0;
my $tsv_type;
my $coreFuncRegion_h;
if (exists $config{CoreFuncRegion}) {
    require CheckInOut;
    $coreFuncRegion_h = CheckInOut->new( db => $config{CoreFuncRegion} );
}

foreach my $var_file (@ARGV) {
    my $varf_h;
    my $tsv_sample_name;
    if ($type ne "tsv") {
        $varf_h = Vcf->new( file => $var_file, version => '4.1' );
        $varf_h->parse_header();
    }
    else {
        my $base_filename = basename($var_file);
        if ($base_filename =~ /^var\-([^\\\/]+)\.tsv(\.bz2)?$/) {
            $tsv_sample_name = $1;
            if (!defined $2) {
                open (CG, $var_file) or confess "Error: [$var_file] $!";
            }
            else {
                open (CG, "bzip2 -dc $var_file |") or confess "Error: [$var_file] $!";
            }
            $varf_h = \*CG;
            $tsv_type = "var";
        }
        elsif ($base_filename =~ /^masterVarBeta\-([^\\\/]+)\.tsv(\.bz2)?$/) {
            $tsv_sample_name = $1;
            if (!defined $2) {
                open (CG, $var_file) or confess "Error: [$var_file] $!";
            }
            else {
                open (CG, "bzip2 -dc $var_file |") or confess "Error: [$var_file] $!";
            }
            $varf_h = \*CG;
            $tsv_type = "master";
        }
        else {
            confess ("Error: cannot recognized the file name of tsv result.\n",
            "      [$var_file] should be in format of var-[ASM_ID].tsv.bz2,\n",
            "      of in format of masterVarBeta-[ASM_ID].tsv.bz2.\n",
            "      ASM_ID is the sample ID in CG's system. which defined in\n",
            "      Standard Sequencing Service Data Format v2.4");
        }
    }

    my $read_stat = 1;
    my $pre_chr;
    my $sample_ids;
    while ($read_stat) {
        my @tobe_anno_vars = ();

        my ($t0, $t1);

        if ( defined $timer ) {
            $t0 = [ gettimeofday ];
        }

        if ( $type ne "tsv" ) {
            ( $read_stat, @tobe_anno_vars ) = read_vcf($varf_h);
        }
        else {
            if ( $tsv_type eq "var" ) {
                ( $read_stat, @tobe_anno_vars ) =
                  read_tsv( $varf_h, $tsv_sample_name );
            }
            elsif ( $tsv_type eq "master" ) {
                ( $read_stat, @tobe_anno_vars ) = 
                  read_master( $varf_h, $tsv_sample_name );
            }
        }

        if ( defined $timer ) {
            $t1 = [ gettimeofday ];
            $read_record_time += tv_interval( $t0, $t1 );
            if ( $read_stat ) {
                $main_record_count ++;
            }
            if ( $main_record_count > $BUFFER_SIZE ) {
                print STDERR "[TIMER RECORD_READ $BUFFER_SIZE] : $read_record_time\n";
                $read_record_time = 0;
                $main_record_count = 0;
            }
        }

        if ( !$read_stat
            or ( defined $pre_chr and $pre_chr ne $tobe_anno_vars[0]->{chr} ) )
        {
            # check homo ref vars to be add before chr changing or end of file
            if (
                    defined $HomoRef_Var
                and defined $sample_ids
                and defined $pre_chr 
                and exists $HomoRef_Var->{$pre_chr}
                and (  !exists $HRVanno_opt->{$pre_chr}
                    or !exists $HRVanno_opt->{$pre_chr}->{$var_file} )
              )
            {
                my @chrom_final_add = ();
                foreach my $var_in_homoR ( @{ $HomoRef_Var->{$pre_chr} } ) {
                    if (   !exists $var_in_homoR->{varfile}
                        or !exists $var_in_homoR->{varfile}->{$var_file} )
                    {
                        if ( my $new_gen =
                            gen_homoRef( $var_in_homoR, $sample_ids ) )
                        {
                            push( @chrom_final_add, $new_gen );
                        }
                        $var_in_homoR->{varfile}->{$var_file} = 1;
                    }
                }

                $enqueue_count += scalar(@chrom_final_add);
                $anno_q->enqueue(@chrom_final_add);
                $HRVanno_opt->{$pre_chr}->{$var_file} = 1;
            }
        }

        
        if (   $enqueue_count >= $BUFFER_SIZE
            or !$read_stat
            or ( defined $pre_chr and $pre_chr ne $tobe_anno_vars[0]->{chr} ) )
        {
            lock($cached_count);

            until ( $cached_count >= $enqueue_count ) {
                my $abs_timeout = time() + $SINGLE_TIMEOUT;
                unless ( cond_timedwait( $cached_count, $abs_timeout ) ) {
                    confess(
                        "Error: TimeOut [$SINGLE_TIMEOUT].",
                        " Please check whether the database index",
                        " works well or maybe all threads have exited."
                    );
                }
            }

            print_cache();

            $enqueue_count -= $cached_count;
            @anno_cache   = ();
            $cached_count = 0;
        }

        if ($read_stat) {

            $total_record_count ++;

            # intilize sample ids in current vcf file
            if ( !defined $sample_ids ) {
                $sample_ids = [ sort keys %{ $tobe_anno_vars[0]->{sample} } ];
            }

            # add reference allele to vcf no-called region due to
            # by default only mutant site will be output in VCF
            check_to_add_homoRef( \@tobe_anno_vars, $var_file, $sample_ids )
              if ( defined $HomoRef_Var );

            $enqueue_count += scalar(@tobe_anno_vars);
            $pre_chr = $tobe_anno_vars[0]->{chr};
            $anno_q->enqueue(@tobe_anno_vars);

            my $qcusing_chr = 'chr'.$pre_chr if ($pre_chr !~ /^chr/);
            if ( defined $in_msqc and defined $msqc_rst and exists $QCsites->{$qcusing_chr} ) {
                # due to the co-presence of variants for one locus
                # we check ms-qc sites here
                foreach my $qcidx (keys %{$$QCsites{$qcusing_chr}}) {
                    my $cur_qc_pos = $$QCsites{$qcusing_chr}{$qcidx}{"pos"};
                    my $cur_qc_ref = $$QCsites{$qcusing_chr}{$qcidx}{"ref"};
                    foreach my $tba (@tobe_anno_vars) {
                        if ( exists $tba->{end} ) {    # CG var format
                            if (    $tba->{begin} < $cur_qc_pos
                                and $tba->{end} >= $cur_qc_pos )
                            {
                                if ( $tba->{variantSequence} eq
                                    $tba->{referenceSequence} )
                                {    # reference called allele
                                    foreach my $samp ( keys %{ $tba->{sample} } ) {
                                        $ngsQCrst->{$samp}->[$qcidx] .= $cur_qc_ref
                                          if ( !exists $ngsQCrst->{$samp}
                                            or !defined $ngsQCrst->{$samp}->[$qcidx]
                                            or $ngsQCrst->{$samp}->[$qcidx] ne
                                            '.' );
                                    }
                                }
                                elsif ($tba->{end} - $tba->{begin} == 1
                                        and $tba->{variantSequence} =~ /^[ACGT]$/i )
                                {
                                    foreach my $samp ( keys %{ $tba->{sample} } ) {
                                        $ngsQCrst->{$samp}->[$qcidx] .=
                                          $tba->{variantSequence}
                                          if ( !exists $ngsQCrst->{$samp}
                                            or !defined $ngsQCrst->{$samp}->[$qcidx]
                                            or $ngsQCrst->{$samp}->[$qcidx] ne
                                            '.' );
                                    }
                                }
                                else {
                                    # other cases all equal to not available
                                    # This will overide any called allele
                                    foreach my $samp ( keys %{ $tba->{sample} } ) {
                                        $ngsQCrst->{$samp}->[$qcidx] = '.';
                                    }
                                }
                            }
                        }
                        else {                         # vcf var format
                            if (
                                $tba->{begin} <= $cur_qc_pos
                                and ( $tba->{begin} +
                                    length( $tba->{referenceSequence} ) -
                                    1 ) >= $cur_qc_pos
                              )
                            {
                                if ( $tba->{variantSequence} eq
                                    $tba->{referenceSequence} )
                                {    # reference called allele
                                    foreach my $samp ( keys %{ $tba->{sample} } ) {
                                        $ngsQCrst->{$samp}->[$qcidx] .= $cur_qc_ref
                                          if ( !exists $ngsQCrst->{$samp}
                                            or !defined $ngsQCrst->{$samp}->[$qcidx]
                                            or $ngsQCrst->{$samp}->[$qcidx] ne
                                            '.' );
                                    }
                                }
                                elsif ( $tba->{begin} == $cur_qc_pos
                                    and $tba->{variantSequence} =~ /^[ACGTN]+$/ )
                                {
                                    # assume anchor sites to be available
                                    # for qc using
                                    my $called_allele =
                                      substr( $tba->{variantSequence}, 0, 1 );
                                    foreach my $samp ( keys %{ $tba->{sample} } ) {
                                        if ( $called_allele ne 'N' ) {
                                            $ngsQCrst->{$samp}->[$qcidx] .=
                                              $called_allele
                                              if ( !exists $ngsQCrst->{$samp}
                                                or !
                                                defined $ngsQCrst->{$samp}->[$qcidx]
                                                or $ngsQCrst->{$samp}->[$qcidx] ne
                                                '.' );
                                        }
                                        else {
                                            $ngsQCrst->{$samp}->[$qcidx] = '.';
                                        }
                                    }
                                }
                                else {    # any ranged or complex variants
                                    foreach my $samp ( keys %{ $tba->{sample} } ) {
                                        $ngsQCrst->{$samp}->[$qcidx] = '.';
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        elsif (defined $in_msqc and defined $msqc_rst) { # finished reading
            # complete qc-sites info for each sample
            foreach my $samp (@$sample_ids) {
                for (my $i = 0; $i < $nQCs; $i++) {
                    foreach my $ch ( sort keys %$QCsites ) {
                        next if (!exists $QCsites->{$ch}->{$i});
                        if (   !exists $ngsQCrst->{$samp}
                            or !defined $ngsQCrst->{$samp}->[$i] )
                        {
                            $ngsQCrst->{$samp}->[$i] = $QCsites->{$ch}->{$i}->{ref};
                        }
                    }
                }
            }
        }
    }

    if ($type ne "tsv") {
        $varf_h->close();
    }
    else {
        close $varf_h;
    }
}
close $OutFp;

# Output all samples qc comparison into a single file
if (defined $msqc_rst and defined $in_msqc) {
    print QCSTAT join("\t", "#SampleID", "Class", 1 .. 21)."\n";
    foreach my $samp (sort keys %$msQCrst) {

        if (!exists $ngsQCrst->{$samp}) {
            print STDERR "[Warning] the MS QC result contains some samples not in input. ($samp) \n";
            next;
        }

        my @qc_cmp_stat = ();
        for (my $i = 0; $i < $nQCs; $i++) {
            if (   $msQCrst->{$samp}->[$i] eq '.'
                or $ngsQCrst->{$samp}->[$i] eq '.'
                or $msQCrst->{$samp}->[$i] =~ /[^ACGT]/
                or $ngsQCrst->{$samp}->[$i] =~ /[^ACGT]/ )
            {
                $qc_cmp_stat[$i] = '?';
            }
            else {
                my $uniform_ms =
                  uc(
                    join( "", sort ( split( //, $msQCrst->{$samp}->[$i] ) ) ) );
                my $uniform_ngs =
                  uc(
                    join( "", sort ( split( //, $ngsQCrst->{$samp}->[$i] ) ) )
                  );
                if ( $uniform_ms eq $uniform_ngs ) {
                    $qc_cmp_stat[$i] = 'ok';
                }
                else {
                    $qc_cmp_stat[$i] = 'unmatch';
                }
            }
        }

        print QCSTAT join("\t", $samp, "(MS)",  @{$msQCrst->{$samp}})."\n";
        print QCSTAT join("\t", $samp, "(NGS)", @{$ngsQCrst->{$samp}})."\n";
        print QCSTAT join("\t", $samp, "(RST)", @qc_cmp_stat)."\n";
    }
    close QCSTAT;
}

$anno_q->insert( 0, ((undef) x ($NUM_THR - 1)) );
sleep 1;
foreach my $thr (threads->list(threads::joinable)) {
    $thr->join();
}
sleep 1;
my @act_thrs = threads->list(threads::running);
foreach my $thr (@act_thrs) {
    $thr->set_thread_exit_only(1);
    $thr->kill('KILL')->detach();
}

# Output Statistic Information
print STDERR "\n\n============= Variation Statistics ============\n";
foreach my $sampleID (sort keys %var_count) {
    my $sample_total = $var_count{$sampleID}{total};
    my $no_call_total = (exists $no_call_var_count{$sampleID}) ? $no_call_var_count{$sampleID} : 0;
    my $real_total = $sample_total + $no_call_total;
    print STDERR "[INFO] $sampleID - called ratio : "
      . "$sample_total / $real_total ("
      . sprintf( "%.3f", ( $sample_total / $real_total * 100 ) ) . "%)\n";

    delete $var_count{$sampleID}{total};
    foreach my $classID ( sort keys %{ $var_count{$sampleID} } ) {
        foreach my $sub_classID ( sort keys %{ $var_count{$sampleID}{$classID} } )
        {
            print STDERR "[INFO] $sampleID - $classID - $sub_classID : "
              . "$var_count{$sampleID}{$classID}{$sub_classID} / $sample_total ("
              . sprintf(
                "%.3f",
                (
                    $var_count{$sampleID}{$classID}{$sub_classID} /
                      $sample_total * 100
                )
              ) . "%)\n";
        }
    }

    $tv{$sampleID} = 0 if ( !exists $tv{$sampleID} );
    $ti{$sampleID} = 0 if ( !exists $ti{$sampleID} );
    print STDERR "[INFO] $sampleID - ti/tv : "
      . "$ti{$sampleID} /  $tv{$sampleID} = "
      . ( ( $tv{$sampleID} == 0 )
        ? "inf"
        : sprintf( "%.3f", $ti{$sampleID} / $tv{$sampleID} ) )
      . "\n";
    if ( exists $sex_varcount{$sampleID} ) {
        $sex_varcount{$sampleID}{het} = 0
          if ( !exists $sex_varcount{$sampleID}{het} );
        print STDERR "[INFO] $sampleID - het in sex : "
          . "$sex_varcount{$sampleID}{het} / $sex_varcount{$sampleID}{total} ("
          . sprintf( "%.3f",
            $sex_varcount{$sampleID}{het} /
              $sex_varcount{$sampleID}{total} *
              100 )
          . "%)\n";
    }

    my $anno_total = $annotation_count{$sampleID}{total};
    delete $annotation_count{$sampleID}{total};
    foreach my $classID ( sort keys %{ $annotation_count{$sampleID} } ) {
        foreach my $sub_classID ( sort keys %{ $annotation_count{$sampleID}{$classID} } )
        {
            print STDERR "[INFO] $sampleID - $classID - $sub_classID : "
              . "$annotation_count{$sampleID}{$classID}{$sub_classID} / $anno_total ("
              . sprintf(
                "%.3f",
                (
                    $annotation_count{$sampleID}{$classID}{$sub_classID} /
                      $anno_total * 100
                )
              ) . "%)\n";
        }
    }
}
print STDERR "== Annotation Done at ". localtime() . " ==\n";

$finished = 1;
sleep 1;
exit 0;


END {
    if (!defined $finished and defined $anno_q) {
        print STDERR "Exception caught, waiting for cleaning threads ";
        $anno_q->insert( 0, ((undef) x ($NUM_THR - 1)) );
        sleep 1;
        foreach my $thr (threads->list(threads::joinable)) {
            $thr->join();
        }
        sleep 1;
        foreach my $thr (threads->list(threads::running)) {
            $thr->set_thread_exit_only(1);
            $thr->exit(1) if $thr->can('exit');
        }
        sleep 1;
        print STDERR "\nAbnormal Returned\n";
        exit 1;
    }
}

# Interrupt Handler
sub quit {
    my $signal = shift;
    print STDERR "Caught signal $signal, waiting for cleaning threads ...";
    if (defined $anno_q) {
        $anno_q->insert( 0, ((undef) x ($NUM_THR - 1)) );
        sleep 1;
        foreach my $thr (threads->list(threads::joinable)) {
            $thr->join();
        }
        sleep 1;
        foreach my $thr (threads->list(threads::running)) {
            $thr->set_thread_exit_only(1);
            $thr->exit(1) if $thr->can('exit');
        }
    }
    print STDERR "\nInterrupted!\n";
    $finished = 1;
    exit 1;
}

sub get_ref {
    my ( $fai_h, $chr, $start, $end ) = @_;
    return "" if ($start >= $end);
    return "=" if ($MAX_REF_LEN > 0 and $end - $start > $MAX_REF_LEN);
    $chr = 'chr' . $chr if ( $chr !~ /^chr/ );
    if ( $chr =~ /^chrM/i ) {
        $chr = 'chrM_NC_012920.1';
    }
    my $region = $chr.':'.($start+1).'-'.$end;
    my $subseq = $fai_h->getseq($region);
    if (!defined $subseq or $subseq eq "") {
        carp "[Warning] Fail to get reference seq for $region\n" if (!defined $quiet);
        return "=";
    }
    return uc($subseq);
}


sub anno_threads {

    my %common_opts = ();
    if (defined $quiet) {
        $common_opts{quiet} = 1;
    }

    require GetHGMD;
    require GetCGDinfo;
    require GetOmim;
    require Faidx;
    require GetClinVar;

    my ($beda, $panel_h, $omim_h, $hgmd_h, $clinVar_h, $vusVar_h, $cgd_h, $hg19_fai);

    $SIG{KILL} = sub {
        $beda->DESTROY() if (defined $beda and $beda->can('DESTROY'));
        $cgd_h->DESTROY() if (defined $cgd_h and $cgd_h->can('DESTROY'));
        $hgmd_h->DESTROY() if (defined $hgmd_h and $hgmd_h->can('DESTROY'));
        $clinVar_h->DESTROY() if (defined $clinVar_h and $clinVar_h->can('DESTROY'));
        $vusVar_h->DESTROY() if (defined $vusVar_h and $vusVar_h->can('DESTROY'));
        $panel_h->DESTROY() if (defined $panel_h and $panel_h->can('DESTROY'));
        $hg19_fai->DESTROY() if (defined $hg19_fai and $hg19_fai->can('DESTROY'));
        carp "Thread ", threads->tid(), ": killed"  if (defined $debug);
        { no warnings 'threads'; cond_signal($cached_count); }
        threads->exit(1);
    };
    local $SIG{__DIE__} = sub {
        my $dead_msg = shift;
        $beda->DESTROY() if (defined $beda and $beda->can('DESTROY'));
        $cgd_h->DESTROY() if (defined $cgd_h and $cgd_h->can('DESTROY'));
        $hgmd_h->DESTROY() if (defined $hgmd_h and $hgmd_h->can('DESTROY'));
        $clinVar_h->DESTROY() if (defined $clinVar_h and $clinVar_h->can('DESTROY'));
        $vusVar_h->DESTROY() if (defined $vusVar_h and $vusVar_h->can('DESTROY'));
        $panel_h->DESTROY() if (defined $panel_h and $panel_h->can('DESTROY'));
        $hg19_fai->DESTROY() if (defined $hg19_fai and $hg19_fai->can('DESTROY'));
        carp "Unexpected Error Meet in Thread ", threads->tid(), ": $dead_msg";
        { no warnings 'threads'; cond_signal($cached_count); }
        threads->exit(1);
    };


    my ($t0, $t1);
    my %thread_timer = ();

    # involke anno engine
    $hg19_fai = Faidx->new( $hg19_fa );
    $beda = BedAnno->new( db => $db, tr => $trDB, %opts, %common_opts );

    if ($config{Panel_Control} ne ".") {
        require GetVcfAF;
        $panel_h = GetVcfAF->new(
            db => $config{Panel_Control}, %common_opts
        );
    }
    
    $omim_h = GetOmim->new( acdb => $omim_data );
    $hgmd_h = GetHGMD->new( db => $hgmd, %common_opts );
    $clinVar_h = GetClinVar->new( db => $clinVar );
    $vusVar_h = GetClinVar->new( db => $vusVar );
    $cgd_h = GetCGDinfo->new( db => $cgd_data, quiet => 1 );

    if ( $write_version and 1 eq threads->tid() ) { # to only execute once
        open (VERLIST, ">", $verList) or die "Error: [$verList] $!";
        foreach my $db (sort keys %DB_VERSION) {
            printf VERLIST "%20s : %s\n", $db, $DB_VERSION{$db};
        }
        close VERLIST;
    }

    my $cur_thread_varcount = 0;
    while (my $rquery_var = $anno_q->dequeue()) {
        my $anno;
        my $reformat_hash;
        print STDERR "[Info from Thread ".threads->tid()."] current var: ".Dumper($rquery_var) if (defined $debug);
        eval {
            local $SIG{__DIE__} = sub {
                die "[Annotation Failed] "
                  . join( "",
                    "[Thread ",
                    threads->tid(),
                    " Var ",
                    "$rquery_var->{chr},",
                    "$rquery_var->{begin},",
                    "$rquery_var->{referenceSequence},",
                    "$rquery_var->{variantSequence}] " )
                  . shift;
            };

            if ($rquery_var->{referenceSequence} eq '=') { # only for tsv format
                $rquery_var->{referenceSequence} = get_ref(
                    $hg19_fai, @$rquery_var{qw{chr begin end}}
                );
            }

            $anno = $beda->anno( $rquery_var );
        };

        if ($@ =~ /Annotation Failed/) {
            print STDERR $@;
            $reformat_hash = convert_null_var($rquery_var);
        }
        else {
            # remove object properties for shared_clone
            $reformat_hash->{var}    = $anno->{var}->TO_JSON();
            
            my $query_var_hash = {
                chr   => $reformat_hash->{var}->{chr},
                start => $reformat_hash->{var}->{pos},
                end   => $reformat_hash->{var}->{end},
                ref   => $reformat_hash->{var}->{ref},
                alt   => $reformat_hash->{var}->{alt},
            };

            if (defined $panel_h and $reformat_hash->{var}->{guess} ne 'no-call') {
                eval {
                    local $SIG{__DIE__} = sub {
                        die "[PanelSQL failed] " 
                          . join( "",
                            "[Thread ",
                            threads->tid(),
                            " Var ",
                            "$query_var_hash->{chr},",
                            "$query_var_hash->{start},",
                            "$query_var_hash->{end},",
                            "$query_var_hash->{ref},",
                            "$query_var_hash->{alt}] " )
                          . shift;
                    };
                    $reformat_hash->{var}->{panelctrl} = $panel_h->getAF(
                        @$query_var_hash{qw(chr start end ref alt)}
                    );
                };
                if ($@ =~ /PanelSQL failed/) {
                    print STDERR $@;
                    $reformat_hash->{var}->{PanelSQL_Fail} = 1;
                }
            }

            eval {
                local $SIG{__DIE__} = sub {
                    die "[clinVar failed] "
                      . join( "",
                        "[Thread ",
                        threads->tid(),
                        " Var ",
                        "$query_var_hash->{chr},",
                        "$query_var_hash->{start},",
                        "$query_var_hash->{end},",
                        "$query_var_hash->{ref},",
                        "$query_var_hash->{alt}] " )
                      . shift;
                };
                my $clinVar_hit = $clinVar_h->getCL(
                    @$query_var_hash{qw(chr start end ref alt)}
                );
                if (%$clinVar_hit) { # hit record in clinVar db
                    $reformat_hash->{var}->{clinVar} = $clinVar_hit;
                }
                else {
                    $reformat_hash->{var}->{clinVar} = $vusVar_h->getCL(
                        @$query_var_hash{qw(chr start end ref alt)}
                    );
                }
            };
            if ( $@ =~ /^clinVar failed/ ) {
                print STDERR $@;
                $reformat_hash->{var}->{clinVar_Fail} = 1;
            }

            eval {
                local $SIG{__DIE__} = sub {
                    die "[get_flanks] "
                      . join( "",
                        "[Thread ",
                        threads->tid(),
                        " Var ",
                        "$query_var_hash->{chr},",
                        "$query_var_hash->{start},",
                        "$query_var_hash->{end},",
                        "$query_var_hash->{ref},",
                        "$query_var_hash->{alt}] " )
                      . shift;
                };
                $reformat_hash->{var}->{flanks} = get_flanks($hg19_fai, $query_var_hash);
            };
            if ( $@ =~ /^\[get_flanks/ ) {
                print STDERR $@;
                $reformat_hash->{var}->{flanks} = ".";
            }

            if ( exists $anno->{trInfo} ) {
                my $major_var_info = $anno->{var}->{varName};
                $reformat_hash->{trInfo} = $anno->{trInfo};
                foreach my $tr ( sort keys %{ $reformat_hash->{trInfo} } ) {

                    my $curtr_ent = $reformat_hash->{trInfo}->{$tr};
                    $curtr_ent->{major} =
                      ( $major_var_info =~ /^$tr/ ) ? "Y" : "N";

                    eval {
                        local $SIG{__DIE__} = sub {
                            die "[GeneSql Fail] "
                              . join( "",
                                "[Thread ",
                                threads->tid(),
                                " for $tr] " )
                              . shift;
                        };

                        $curtr_ent->{omim} =
                          $omim_h->getAnnoComb( $curtr_ent->{geneId} );
                        $curtr_ent->{cgd} =
                          $cgd_h->getCGD( $curtr_ent->{geneId} );
                    };
                    if ( $@ =~ /GeneSql Fail/ ) {
                        print STDERR $@;
                        $curtr_ent->{omim} = {};
                        $curtr_ent->{cgd}  = {};
                        $curtr_ent->{GeneSQL_Fail} = 1;
                    }

                    $curtr_ent->{TCD} =
                      ( defined $GeneTestCode and exists $GeneTestCode->{ $curtr_ent->{geneSym} } )
                      ? $GeneTestCode->{ $curtr_ent->{geneSym} }
                      : ".";

                    # There's no need to query tr-extra database for no-call var
                    next if ($reformat_hash->{var}->{guess} eq 'no-call');

                    my $trQuery = {%$query_var_hash};

                    @$trQuery{qw(
                        transcript geneid genesym cHGVS 
                    )} = ($tr, @$curtr_ent{qw(
                        geneId geneSym c 
                    )});

                    $trQuery->{pHGVS} = $curtr_ent->{p} if (exists $curtr_ent->{p} and defined $curtr_ent->{p});
                    $trQuery->{pHGVS3} = $curtr_ent->{p3} if (exists $curtr_ent->{p3} and defined $curtr_ent->{p3});

                    my $query_alt_opt = 0;
                    my $altTrQ = {%$trQuery};
                    if (exists $curtr_ent->{alt_cHGVS} or exists $curtr_ent->{alt_pHGVS}) {
                        $query_alt_opt = 1;
                        $altTrQ->{cHGVS} = $curtr_ent->{alt_cHGVS} if (exists $curtr_ent->{alt_cHGVS});
                        $altTrQ->{pHGVS} = $curtr_ent->{alt_pHGVS} if (exists $curtr_ent->{alt_pHGVS});
                        $altTrQ->{pHGVS3} = $curtr_ent->{alt_p3} if (exists $curtr_ent->{alt_p3});
                    }

                    my $query_std_opt = 0;
                    my $stdTrQ = {%$trQuery};
                    if (exists $curtr_ent->{standard_cHGVS} or exists $curtr_ent->{standard_pHGVS}) {
                        $query_std_opt = 1;
                        $stdTrQ->{cHGVS} = $curtr_ent->{standard_cHGVS} if (exists $curtr_ent->{standard_cHGVS});
                        $stdTrQ->{pHGVS} = $curtr_ent->{standard_pHGVS} if (exists $curtr_ent->{standard_pHGVS});
                        $stdTrQ->{pHGVS3} = $curtr_ent->{standard_p3} if (exists $curtr_ent->{standard_p3});
                    }

                    if (    exists $curtr_ent->{prRef}
                        and exists $curtr_ent->{prAlt}
                        and length( $curtr_ent->{prRef} ) == 1
                        and length( $curtr_ent->{prAlt} ) == 1 )
                    {
                        $trQuery->{aaref} = $curtr_ent->{prRef};
                        $trQuery->{aaalt} = $curtr_ent->{prAlt};
                    }

                    eval {
                        local $SIG{__DIE__} = sub {
                            die "[localhgmd Fail] "
                              . join( "",
                                "[Thread ",
                                threads->tid(),
                                " $tr: ",
                                $trQuery->{cHGVS},
                                " ($trQuery->{pHGVS})] " )
                              . shift;
                        };
                        $curtr_ent->{hgmd} = $hgmd_h->getHGMD($trQuery);

                        if ( $query_std_opt and 0 >= scalar @{$curtr_ent->{hgmd}} ) {
                            $curtr_ent->{hgmd} = $hgmd_h->getHGMD($stdTrQ);
                        }
                        if ( $query_alt_opt and 0 >= scalar @{$curtr_ent->{hgmd}} ) {
                            $curtr_ent->{hgmd} = $hgmd_h->getHGMD($altTrQ);
                        }
                    };
                    if ($@ =~ /localhgmd Fail/) {
                        print STDERR $@;
                        $curtr_ent->{hgmd} = [];
                        $curtr_ent->{LocalHGMD_Fail} = 1;
                    }
                }
            }
        }
        
        $cur_thread_varcount ++;

        print STDERR "[Info from Thread ".threads->tid()."] var after anno:".Dumper($reformat_hash) if (defined $debug);
        {
            lock($cached_count);
            push (@anno_cache, shared_clone({ %$reformat_hash }));
            $cached_count ++;
            cond_signal($cached_count);
        }

        if ( $cur_thread_varcount >= $BUFFER_SIZE ) {
            if ( defined $timer ) {
                foreach my $timekey (sort keys %thread_timer) {
                    print STDERR "[TIMER ".threads->tid()."] $timekey : $thread_timer{$timekey}\n";
                    $thread_timer{$timekey} = 0;
                }
            }
            $cur_thread_varcount = 0;
        }
    }

    $beda->DESTROY() if (defined $beda and $beda->can('DESTROY'));
    $cgd_h->DESTROY() if (defined $cgd_h and $cgd_h->can('DESTROY'));
    $hgmd_h->DESTROY() if (defined $hgmd_h and $hgmd_h->can('DESTROY'));
    $vusVar_h->DESTROY() if (defined $vusVar_h and $vusVar_h->can('DESTROY'));
    $clinVar_h->DESTROY() if (defined $clinVar_h and $clinVar_h->can('DESTROY'));
    $hg19_fai->DESTROY() if (defined $hg19_fai and $hg19_fai->can('DESTROY'));
    return;
}

sub get_flanks {
    my $fa_fai = shift;
    my $qvar = shift;
    my $chr = $qvar->{chr};
    $chr = 'chr'.$chr if ($chr !~ /^chr/);
    if ($chr eq 'chrMT') {
        $chr = 'chrM_NC_012920.1';
    }
    my $left_flank_region =
        $chr . ':'
      . ( $qvar->{start} - $config{FLANK_LEN} + 1 ) . '-'
      . $qvar->{start};
    my $right_flank_region =
        $chr . ':'
      . ( $qvar->{end} + 1 ) . '-'
      . ( $qvar->{end} + $config{FLANK_LEN} );
    my @lr_flanks = ();
    foreach my $rgn ($left_flank_region, $right_flank_region) {
        my $cut_off_seq = $fa_fai->getseq($rgn);
        if (!defined $cut_off_seq or $cut_off_seq eq "") {
            carp "Warning: [$rgn] not available." if (!defined $quiet);
            push (@lr_flanks, "[NA]");
        }
        else {
            push (@lr_flanks, $fa_fai->getseq($rgn));
        }
    }
    my $flankseq = join(".", @lr_flanks);
    return $flankseq;
}

sub convert_null_var {
    my $rquery_var = shift;

    my $anno_error_hash;
    $anno_error_hash->{var} = $rquery_var;

    $anno_error_hash->{var}->{ref} = $rquery_var->{referenceSequence};
    $anno_error_hash->{var}->{alt} = $rquery_var->{variantSequence};
    delete $anno_error_hash->{var}->{referenceSequence};
    delete $anno_error_hash->{var}->{variantSequence};

    $anno_error_hash->{var}->{guess} = "annotation-error";

    if ( exists $rquery_var->{end} ) {    # CG
        $anno_error_hash->{var}->{pos} = $rquery_var->{begin};
    }
    else {
        $anno_error_hash->{var}->{pos} = $rquery_var->{begin} - 1;
        $anno_error_hash->{var}->{end} = $anno_error_hash->{var}->{pos} +
          length( $anno_error_hash->{var}->{ref} );
    }
    delete $anno_error_hash->{var}->{begin};
    return $anno_error_hash;
}

sub judge_titv {
    my ($r, $a) = @_;
    my $to_check = join( "", sort ( $r, $a ) );
    if ($to_check eq 'AG' or $to_check eq 'CT') {
        return 1; # Ti
    }
    else {
        return 0; # Tv
    }
}

# current give filter assessment by
# judging the AD and PLtag, besides possible NB for VCF
# for tsv format(CG) using varFilter column
sub give_filter {
    my $sample_h_in_var = shift;
    if ($type eq "tsv") {
        # CG data
        if (exists $sample_h_in_var->{VF}) {
            if ($sample_h_in_var->{VF} =~ /VQLOW/i) {
                $sample_h_in_var->{filterTag} = 'FAIL';
            }
            elsif ($sample_h_in_var->{VF} =~ /AMBIGUOUS/i) {
                $sample_h_in_var->{filterTag} = 'DUBIOUS';
            }
            elsif ($sample_h_in_var->{VF} =~ /^PASS$/i) {
                $sample_h_in_var->{filterTag} = 'PASS';
            }
            else {
                $sample_h_in_var->{filterTag} = '.';
            }
        }
        else {
            $sample_h_in_var->{filterTag} = '.';
        }
    }
    else {
        # VCF data
        if (   !exists $sample_h_in_var->{AD}
            or $sample_h_in_var->{AD} eq '.'
            or !exists $sample_h_in_var->{PLtag}
            or $sample_h_in_var->{PLtag} eq '.' )
        {
            $sample_h_in_var->{filterTag} = '.';
        }
        else {
            if (   $sample_h_in_var->{AD} <= $AD_DN_THRESHOLD
                or $sample_h_in_var->{PLtag} < 0 )
            {
                $sample_h_in_var->{filterTag} = 'FAIL';
            }
            elsif (
                    $sample_h_in_var->{AD} >= $AD_UP_THRESHOLD
                and $sample_h_in_var->{PLtag} > 0
                and $sample_h_in_var->{VF} eq '.'
                and ( !exists $sample_h_in_var->{NB}
                    or $sample_h_in_var->{NB} eq '.' )
              )
            {
                $sample_h_in_var->{filterTag} = 'PASS';
            }
            else {
                $sample_h_in_var->{filterTag} = 'DUBIOUS';
            }
        }
    }
    return $sample_h_in_var;
}

# check PL score
# return PLtag:
# -1:   FAIL
#  0:   DUBIOUS
#  1:   PASS
sub check_PL {
    my ($PL, $a1, $a2) = @_;
    my @PLs = split(/,/, $PL);
    my $PLtag = 1;
    my $PLindex = pIndex_cal( sort {$a<=>$b} ( $a1, $a2 ) );
    my $gtype_PL = splice(@PLs, $PLindex, 1);
    if ($gtype_PL eq ".") {
        $PLtag = 1;
    }
    elsif ($gtype_PL > $PL_THRESHOLD) {
        $PLtag = -1;
    }
    elsif ($gtype_PL > 0) {
        $PLtag = 0;
    }
    else {
        foreach my $otherPL (@PLs) {
            if ($otherPL eq ".") {
                next;
            }
            elsif ($otherPL == 0) {
                $PLtag = -1;
                last;
            }
            elsif ($otherPL < $PL_UP_THRESHOLD) {
                $PLtag = 0;
            }
        }
    }
    return $PLtag;
}

# calculate the index in GL/PL string, when geno pair is from $p1 and $p2 alleles
sub pIndex_cal {
    my ($p1, $p2) = @_;
    return ($p2 * ($p2 + 1) / 2 + $p1);
}

# vcf reader, input a vcf handler, return a var array.
# IMPRECISE variant will be ignored.
sub read_vcf {
    my $vcf_h = shift;
    my $rec = $vcf_h->next_data_hash();
    return (0) unless ($rec);
    return read_vcf($vcf_h) if ( exists $$rec{INFO}{IMPRECISE} );

    my $filter_tag = join(";", @{$$rec{FILTER}});
    $filter_tag =~ s/PASS/./ig;

    my ( $chr, $vcf_pos, $vcf_ref, @vcf_alts ) =
      ( $$rec{CHROM}, $$rec{POS}, $$rec{REF}, @{ $$rec{ALT} } );

    $chr =~ s/^chr//i;
    
    # unshift vcf ref to all case
    my @all_case = ($vcf_ref, @vcf_alts);

    my $nocall_var = {
        chr               => $chr,
        begin             => ( $vcf_pos - 1 ),
        end               => ( $vcf_pos + length($vcf_ref) - 1 ),
        referenceSequence => $vcf_ref,
        variantSequence   => '?',
#        _INFO             => $$rec{INFO},
    };

    my @vars = ();
    for my $i ( 0 .. $#all_case ) {
        # ignore the non precise variation for vcf format
        # ignore large sv or breakend for vcf format
        return read_vcf($vcf_h)
          if ( $all_case[$i] =~ /[<>\[\]]/
            or ( $all_case[$i] ne "." and $all_case[$i] =~ /\./ ) );

        $all_case[$i] = $all_case[0] if ($all_case[$i] eq '.');

        # this is just indicate it's from vcf
        # DO NOT add "end" key or you will have to
        # change begin to 0-based position
        my $var = {
            chr               => $chr,
            begin             => $vcf_pos,
            referenceSequence => $vcf_ref,
            variantSequence   => $all_case[$i],
#            _INFO             => $$rec{INFO},
        };

        push (@vars, $var);
    }

    # assign sample info to var
    # "sample" => {
    #   $sample => {
    #       zygosity => $zygosity,
    #       AD => $allelicDepth,
    #       AR => $allelicRatio,
    #       VF => $CallerFilterTag,
    #       AI => $AllelicIndex
    #       PLOIDY => $ploidy,
    #       PLTAG => $judge_reliablity_from_PL,
    #    
    #       # Here is customized keys
    #       PB => $Phased_Group_ID,
    #       NB => $Neighbour_Group_ID,
    #   },
    #   ...
    # }

    my $involve_nocall = 0;
    foreach my $sample (sort keys %{$$rec{gtypes}}) {
        my $sample_hash = $$rec{gtypes}{$sample};
        $sample_hash->{GT} =~ $$vcf_h{regex_gt}
          or confess "Could not parse gtype string [$sample_hash->{GT}] [$chr:$vcf_pos]";
        my ($a1, $sep, $a2) = ($1, $2, $3);
        my ($Zyg, $AD, $AR, $PLtag, $PB, $NB) = ('.') x 6;
        $PB = $sample_hash->{PB}
          if ( exists $sample_hash->{PB} and defined $sample_hash->{PB} );
        $NB = $sample_hash->{NB}
          if ( exists $sample_hash->{NB} and defined $sample_hash->{NB} );
        my $AI = 0; # default homo/hap
        my $PLOIDY = (defined $a2 and $a2 ne '') ? 2 : 1;

        my @ADs = ();
        my $total_Dp = 0;
        if (    exists $sample_hash->{AD}
            and defined $sample_hash->{AD}
            and $sample_hash->{AD} ne '.' )
        {
            # no calling case will non-exists the tag.
            @ADs = split( /,/, $sample_hash->{AD} );
            $total_Dp += $_ foreach (@ADs);
            $pre_total_AD{$sample} = $total_Dp;
        }

        if (exists $sample_hash->{PL} and defined $a2 and $a2 ne '') {
            $PLtag = check_PL($sample_hash->{PL}, $a1, $a2);
        }

        if (!defined $a2 or $a2 eq '' or $a1 eq $a2) {
            $AI = 0;
            my $cur_sample_h;
            if ($a1 eq '.') {
                $nocall_var->{sample}{$sample} = {};
                $cur_sample_h = $nocall_var->{sample}{$sample};
                $involve_nocall = 1;
                $Zyg = 'no-call';
            }
            else {
                $vars[$a1]->{sample}{$sample} = {};
                $cur_sample_h = $vars[$a1]->{sample}{$sample};
                if (defined $a2 and $a2 eq '0') {
                    $Zyg = 'hom-ref';
                }
                elsif (defined $a2 and $a2 ne '') {
                    $Zyg = 'hom-alt';
                }
                elsif ($a1 eq '0') {
                    $Zyg = 'hap-ref';
                }
                else {
                    $Zyg = 'hap-alt';
                }

                if ($total_Dp > 0) {
                    if (!defined $ADs[$a1]) {
                        carp "AD index error! [$chr, $vcf_pos, $sample]"
                    }
                    else {
                        $AD = $ADs[$a1];
                        $AR = sprintf("%.2f", $AD/$total_Dp);
                    }
                }
            }
            $cur_sample_h->{zygosity} = $Zyg;
            $cur_sample_h->{NB} = $NB;
            $cur_sample_h->{PB} = $PB;
            $cur_sample_h->{AD} = $AD;
            $cur_sample_h->{AR} = $AR;
            $cur_sample_h->{VF} = $filter_tag;
            $cur_sample_h->{AI} = $AI;
            $cur_sample_h->{PLOIDY} = $PLOIDY;
            $cur_sample_h->{PLtag} = $PLtag;
            $cur_sample_h = give_filter($cur_sample_h);
        }
        else {
            my @two_sample_h = ();
            my @aids = ($a1, $a2);
            for my $i ( 0, 1 ) {
                ($Zyg, $AD, $AR) = ('.') x 3;
                my $cur_aid  = $aids[$i];
                my $pair_aid = $aids[ $i - 1 ];
                $AI = $i + 1;
                if ($cur_aid eq '.') {
                    $nocall_var->{sample}{$sample} = {};
                    $two_sample_h[$i] = $nocall_var->{sample}{$sample};
                    $involve_nocall = 1;
                    $Zyg = 'no-call';
                }
                else {
                    $vars[$cur_aid]->{sample}{$sample} = {};
                    $two_sample_h[$i] = $vars[$cur_aid]->{sample}{$sample};
                    if ($pair_aid eq '.') {
                        $Zyg = ($cur_aid eq '0') ? 'half-ref' : 'half-alt';
                    }
                    elsif ($pair_aid ne '0') {
                        $Zyg = 'het-alt';
                    }
                    else { # pair_aid eq '0'
                        $Zyg = 'het-ref';
                    }
                    if ($total_Dp > 0) {
                        if (!defined $ADs[$cur_aid]) {
                            carp "AD index error! [$chr, $vcf_pos, $sample]"
                        }
                        else {
                            $AD = $ADs[$cur_aid];
                            $AR = sprintf("%.2f", $AD/$total_Dp);
                        }
                    }
                }
                $two_sample_h[$i]->{zygosity} = $Zyg;
                $two_sample_h[$i]->{NB} = $NB;
                $two_sample_h[$i]->{PB} = $PB;
                $two_sample_h[$i]->{AD} = $AD;
                $two_sample_h[$i]->{AR} = $AR;
                $two_sample_h[$i]->{AI} = $AI;
                $two_sample_h[$i]->{PLOIDY} = $PLOIDY;
                $two_sample_h[$i]->{VF} = $filter_tag;
                $two_sample_h[$i]->{PLtag} = $PLtag;
                $two_sample_h[$i] = give_filter($two_sample_h[$i]);
            }
        }
    }

    if ($involve_nocall) {
        push (@vars, $nocall_var);
    }

    my @available_vars = ();
    my $refVar = 1; 
    foreach my $v (@vars) {
        my $zyg_tmp;
        if ($refVar) {
            $zyg_tmp = 'ref';
            $refVar = 0;
        }
        else {
            $zyg_tmp = 'alt';
        }
        if (!exists $v->{sample}) { # those dabase like vcf don't have sample infos
            $v->{sample} = {
                "nullSample" => {
                    zygosity => 'homo-'.$zyg_tmp,
                    AD => '.',
                    AR => '.',
                    VF => '.',
                    AI => 0,
                    PLOIDY => 2,
                    PLtag => '.',
                    PB => '.',
                    NB => '.',
                    filterTag => '.',
                }
            }
        }
    }
    return (1, @vars);
}

# check to add homo Ref vars to variant list
sub check_to_add_homoRef {
    my ($ready_anno_vars, $cur_file, $cur_samples) = @_;
    my $test_var = $ready_anno_vars->[0];
    my $chr = $test_var->{chr};
    if (
        exists $HomoRef_Var->{$chr}
        and (  !exists $HRVanno_opt->{$chr}
            or !exists $HRVanno_opt->{$chr}->{$cur_file} )
      )
    {
        my $start = $test_var->{begin};
        $start -= 1 if ( $type ne 'tsv' );
        my $ref = $test_var->{referenceSequence};
        my $end = ( $start + length($ref) );

        my $add_all_opt = 1;
        my @add_before = ();
        foreach my $var_in_homoR ( @{ $HomoRef_Var->{$chr} } ) {
            if ( $var_in_homoR->{begin} >= $end ) {
                $add_all_opt = 0;
                last;
            }
            elsif (!exists $var_in_homoR->{varfile}
                or !exists $var_in_homoR->{varfile}->{$cur_file} )
            {
                if ( $var_in_homoR->{end} <= $start ) {

                    if ( my $new_gen =
                        gen_homoRef( $var_in_homoR, $cur_samples ) )
                    {
                        push( @add_before, $new_gen );
                    }
                    $var_in_homoR->{varfile}->{$cur_file} = 1;
                }
                elsif ( $var_in_homoR->{begin} == $start
                    and $var_in_homoR->{end} == $end )
                {
                    # convert no-call to ref for identical match
                    $ready_anno_vars = fix_no_call_to_ref($ready_anno_vars);
                    $var_in_homoR->{varfile}->{$cur_file} = 1;
                }
                elsif ( $var_in_homoR->{begin} < $end ) {

                    # skip overlapped case
                    $var_in_homoR->{varfile}->{$cur_file} = 1;
                }
            }
        }
        unshift (@$ready_anno_vars, @add_before);

        if ( $add_all_opt == 1 ) {
            $HRVanno_opt->{$chr}->{$cur_file} = 1;
        }
    }
    return $ready_anno_vars;
}

sub fix_no_call_to_ref {
    my $rReady_vars = shift;
    my $ref_var_idx = -1;
    my $no_call_var_idx = -1;
    for (my $i = 0; $i < scalar(@$rReady_vars); $i++) {
        my $cur_var = $rReady_vars->[$i];
        if ($cur_var->{variantSequence} eq '?') {
            $no_call_var_idx = $i;
        }
        if ($cur_var->{referenceSequence} eq $cur_var->{variantSequence}) {
            $ref_var_idx = $i;
        }
    }
    if ($no_call_var_idx < 0) {
        return $rReady_vars;
    }
    elsif ($ref_var_idx < 0) { # no ref-call then fix no-call var
        my $no_call_var = $rReady_vars->[$no_call_var_idx];
        $no_call_var->{variantSequence} = $no_call_var->{referenceSequence};
        foreach my $smp (keys %{$no_call_var->{sample}}) {
            $no_call_var->{sample}->{$smp}->{zygosity} = 'hom-ref';
            my $preAD = $pre_total_AD{$smp}
              if ( exists $pre_total_AD{$smp} );
            $preAD ||= 20;
            $no_call_var->{sample}->{$smp}->{AD} = $preAD;
            $no_call_var->{sample}->{$smp}->{AR} = '1.00';
            $no_call_var->{sample}->{$smp}->{PLtag} = 1;
            $no_call_var->{sample}->{$smp}->{filterTag} = 'PASS';
            $no_call_var->{sample}->{$smp}->{AutoCurated} = 1;
        }
        return $rReady_vars;
    }
    else { # with no-call ref-all, correct ref-all and delete-no-call
        my $no_call_var = $rReady_vars->[$no_call_var_idx];
        my $ref_var = $rReady_vars->[$ref_var_idx];
        foreach my $smp (keys %{$no_call_var->{sample}}) {
            $ref_var->{sample}->{$smp} = $no_call_var->{sample}->{$smp};
            $ref_var->{sample}->{$smp}->{zygosity} = 'hom-ref';
            my $preAD = $pre_total_AD{$smp}
              if ( exists $pre_total_AD{$smp} );
            $preAD ||= 20;
            $ref_var->{sample}->{$smp}->{AD} = $preAD;
            $ref_var->{sample}->{$smp}->{AR} = '1.00';
            $ref_var->{sample}->{$smp}->{PLtag} = 1;
            $ref_var->{sample}->{$smp}->{filterTag} = 'PASS';
            $ref_var->{sample}->{$smp}->{AutoCurated} = 1;
        }
        splice( @$rReady_vars, $no_call_var_idx, 1);
        return $rReady_vars;
    }
}

# generate var entry for homo ref variant
# if all samples specified in sample_ids
# have been annotated on the homoRef_var
# then return undef, otherwise generate 
# a var entry ready to push to be annotated
sub gen_homoRef {
    my ($homoRef_var, $sample_ids) = @_;

    my $return_homoRef = {%$homoRef_var};
    delete $return_homoRef->{sample} if (exists $return_homoRef->{sample});

    foreach my $sample (@$sample_ids) {
        next
          if (  exists $homoRef_var->{sample}
            and exists $homoRef_var->{sample}->{$sample} );
        $homoRef_var->{sample}->{$sample} = 1;

        my $preAD = $pre_total_AD{$sample} if (exists $pre_total_AD{$sample});
        $preAD ||= 20;
        $return_homoRef->{sample}->{$sample} = {
            zygosity  => 'hom-ref',
            NB        => '.',
            PB        => '.',
            AD        => $preAD,
            AR        => '1.00',
            AI        => 0,
            PLOIDY    => 2,
            VF        => '.',
            PLtag     => 1,
            filterTag => 'PASS',
            AutoCurated => 1,
        };
    }

    if (exists $return_homoRef->{sample}) {
        return $return_homoRef;
    }
    else {
        return undef;
    }
}

# tsv file should be carefully sorted
# and formatted. 
# Here is for masterVarBeta-[ASM_ID].tsv.bz2
sub read_master {
    my ($fh, $sample) = @_;
    my $line = <$fh>;
    while ( defined $line
        and ( $line =~ /^#/ or $line =~ /^>/ or $line =~ /^\s*$/ ) )
    {
        $line = <$fh>;
    }
    return (0) if (!defined $line);

    # headers
    # 1. locus
    # 2. ploidy
    # 3. chromosome
    # 4. begin
    # 5. end
    # 6. zygosity
    # 7. varType
    # 8. reference
    # 9. allele1Seq
    # 10. allele2Seq
    # 11. allele1VarScoreVAF
    # 12. allele2VarScoreVAF
    # 13. allele1VarScoreEAF
    # 14. allele2VarScoreEAF
    # 15. allele1VarFilter
    # 16. allele2VarFilter
    # 17. allele1HapLink
    # 18. allele2HapLink
    # 19. allele1XRef
    # 20. allele2XRef
    # 21. allele1Freq
    # 22. allele2Freq
    # 23. allele1AlternativeCalls
    # 24. allele2AlternativeCalls
    # 25. evidenceIntervalId
    # 26. allele1ReadCount
    # 27. allele2ReadCount
    # 28. referenceAlleleReadCount
    # 29. totalReadCount
    # 30. allele1Gene
    # 31. allele2Gene
    # 32. pfam
    # 33. miRBaseId
    # 34. repeatMasker
    # 35. segDupOverlap
    # 36. relativeCoverageDiploid
    # 37. calledPloidy
    # 38. relativeCoverageNondiploid
    # 39. calledLevel
    # 40. bestLAFsingle
    # 41. lowLAFsingle
    # 42. highLAFsingle

    chomp($line);
    my @itms = split(/\t/, $line, 50);

    if ($itms[7] eq '=' and $itms[8] eq '=' and $itms[5] =~ /hom|hap/) {
        # skip this kind of reference call
        # due to not reparse reference bases from fasta
        return read_master( $fh, $sample );
    }
    
    # uniform record
    $itms[2] =~ s/^chr//i;
    my $ploidy = $itms[1];
    my @vars = ();

    if ($itms[8] ne '=') {
        my %a1_info = ();
        @a1_info{qw(
            LocID VAF EAF AltCall refAD totalAD 
        )} = @itms[0, 10, 12, 22, 27, 28 ];

        my $var1 = {
            chr               => $itms[2],
            begin             => $itms[3],
            end               => $itms[4],
            referenceSequence => $itms[7],
            variantSequence   => $itms[8],
    #       _CGMISC           => \%a1_info,
        };


        my $ai_1 = ($itms[5] eq 'no-call' or $ploidy eq '1' or $itms[8] eq $itms[9]) ? 0 : 1;
        my $pb_1 = ($itms[16] =~ /^\d+$/) ? $itms[16] : '.';
        my $ad_1 = ($itms[25] =~ /^\d+$/) ? $itms[25] : 0;
        my $ar_1 =
          ( $a1_info{totalAD} =~ /^\d+$/ and $a1_info{totalAD} > 0 )
          ? sprintf( "%.3f", $ad_1 / $a1_info{totalAD} )
          : '.';
        my $vf_1 = ($itms[14] ne "") ? $itms[14] : '.';

        my $sample_h1 = {
            PB => $pb_1,
            AD => $ad_1,
            AR => $ar_1,
            AI => $ai_1,
            VF => $vf_1,
            PLOIDY => $ploidy,
        };

        my $zy_1 = $itms[5];

        if ($itms[8] eq '?' or $itms[8] eq 'N') {
            $zy_1 = 'no-call';
        }
        elsif ($zy_1 =~ /hap|hom|half/) {
            $zy_1 .= ( $itms[7] eq $itms[8] ) ? '-ref' : '-alt';
        }
        elsif ($zy_1 =~ /het-ref/ and $itms[7] eq $itms[8]) {
            $zy_1 = 'het-alt';
        }
        
        $sample_h1->{zygosity} = $zy_1;
        $sample_h1 = give_filter($sample_h1);
        $var1->{sample}->{$sample} = $sample_h1;
        push (@vars, $var1);
    }

    if ($ploidy eq '2' and $itms[9] ne '=' and $itms[5] !~ /no\-call|hom|hap/) {
        my %a2_info = ();
        @a2_info{qw(
            LocID VAF EAF AltCall refAD totalAD 
        )} = @itms[0, 11, 13, 23, 27, 28];

        my $var2 = {
            chr               => $itms[2],
            begin             => $itms[3],
            end               => $itms[4],
            referenceSequence => $itms[7],
            variantSequence   => $itms[9],
#           _CGMISC           => \%a2_info,
        };

        my $ai_2 = 2;
        my $pb_2 = ($itms[17] =~ /^\d+$/) ? $itms[17] : '.';
        my $ad_2 = ($itms[26] =~ /^\d+$/) ? $itms[26] : 0;
        my $ar_2 =
          ( $a2_info{totalAD} =~ /^\d+$/ and $a2_info{totalAD} > 0 )
          ? sprintf( "%.3f", $ad_2 / $a2_info{totalAD} )
          : '.';
        my $vf_2 = ($itms[15] ne "") ? $itms[15] : '.';

        my $sample_h2 = {
            PB => $pb_2,
            AD => $ad_2,
            AR => $ar_2,
            AI => $ai_2,
            VF => $vf_2,
            PLOIDY => $ploidy,
        };

        my $zy_2;
        if ($itms[9] eq '?' or $itms[9] eq 'N') {
            $zy_2 = 'no-call';
        }
        elsif ($itms[7] eq $itms[8]) {
            $zy_2 = 'het-ref';
        }
        else { # half call case will empty the allele2
            $zy_2 = 'het-alt';
        }
        
        $sample_h2->{zygosity} = $zy_2;
        $sample_h2 = give_filter($sample_h2);
        $var2->{sample}->{$sample} = $sample_h2;
        push (@vars, $var2);
    }

    return (1, @vars);
}

# tsv file should be carefully sorted
# and formatted. 
# Here is for var-[ASM_ID].tsv.bz2
sub read_tsv {
    my ($fh, $sample) = @_;
    my $line = <$fh>;
    while ( defined $line
        and ( $line =~ /^#/ or $line =~ /^>/ or $line =~ /^\s*$/ ) )
    {
        $line = <$fh>;
    }
    return (0) if (!defined $line);

    # headers
    # 1. locus
    # 2. ploidy
    # 3. alleleIndex (all, 1, 2)
    # 4. chr
    # 5. begin
    # 6. end
    # 7. varType
    # 8. ref
    # 9. Call
    # 10. varscoreVAF
    # 11. varscoreEAF
    # 12. varFilter
    # 13. haplink
    # 14. xRef
    # 15. AF
    # 16. altCall

    chomp($line);
    my @itms = split(/\t/, $line, 17);
    my @locus_group = ();
    my $cur_var = tsv_var_parser(\@itms, $sample);
    my $cur_locus = $cur_var->{locus};
    if ( exists $prev_var{$sample} ) {
        # not the first var for this sample
        if ( $prev_var{$sample}->{locus} eq $itms[0] ) {
            # same to current reading
            push( @locus_group, { %{ $prev_var{$sample} } } );
        }
        else {
            if ( $prev_var{$sample}->{referenceSequence} eq
                    $prev_var{$sample}->{variantSequence}
                and $prev_var{$sample}->{referenceSequence} eq '=' )
            {
                # current protocol will skip this kind of
                # reference call, due to not reparse reference
                # bases from fasta.
                $prev_var{$sample} = $cur_var;
                return read_tsv( $fh, $sample );
            }
            my $return_var = { %{ $prev_var{$sample} } };
            if ($return_var->{sample}->{$sample}->{AI} eq '0' or $return_var->{sample}->{$sample}->{PLOIDY} == 1) {
                if (   $return_var->{callerVarType} =~ /^no\-/ )
                {
                    $return_var->{sample}->{$sample}->{zygosity} = 'no-call';
                }
                else {
                    my ($pre, $post) = ();
                    if ($return_var->{sample}->{$sample}->{PLOIDY} == 1) {
                        $pre = "hap";
                    }
                    else {
                        $pre = "hom";
                    }

                    if ($return_var->{callerVarType} =~ /ref/) {
                        $post = 'ref';
                    }
                    else {
                        $post = 'alt';
                    }

                    $return_var->{sample}->{$sample}->{zygosity} = $pre.'-'.$post;
                }

                $prev_var{$sample} = $cur_var;
                return (1, $return_var);
            }
            else {
                # not complete locus?
                carp "[Critical Warning] not complete locus $return_var->{locus} for sample $sample. skipping ...\n";
                $prev_var{$sample} = $cur_var;
                return read_tsv( $fh, $sample );
            }
        }
    }

    while ( $cur_locus eq $cur_var->{locus} ) {
        push (@locus_group, $cur_var);
        my $newline = <$fh>;
        last unless ( $newline );
        my @new_itms = split(/\t/, $newline, 17);
        $cur_var = tsv_var_parser( \@new_itms, $sample );
    }

    if ( $cur_locus ne $cur_var->{locus} ) {
        $prev_var{$sample} = $cur_var;
    }
    my $ref_locus_group = locus_group_parser( \@locus_group );
    return (1, @$ref_locus_group);
}

sub locus_group_parser {
    my $rlocus_group = shift;

    # only 1 sample in tsv's var hash
    my ( $sample ) = keys %{$rlocus_group->[0]->{sample}};

    my %locus_strand_index = ();
    for ( my $i = 0; $i < @$rlocus_group; $i++ ) {
        my $allele_index = $rlocus_group->[$i]->{sample}->{$sample}->{AI};
        if (!exists $locus_strand_index{$allele_index}) {
            $locus_strand_index{$allele_index} = [];
        }
        push ( @{$locus_strand_index{$allele_index}}, $i );
    }

    my ($rlg, $locus_strand_index) = check_combin_group( $rlocus_group, \%locus_strand_index );

    foreach my $ai (sort {$a <=> $b} keys %$locus_strand_index) {
        foreach my $var_idx (@{$locus_strand_index->{$ai}}) {
            my $ai_var = $rlg->[$var_idx];
            if (   $ai_var->{variantSequence} eq '?' )
            {
                $ai_var->{sample}->{$sample}->{zygosity} = 'no-call';
            }
            else {
                my $cur_post =
                  ( $ai_var->{callerVarType} =~ /ref/ ) ? 'ref' : 'alt';
                if ($ai eq '0') {
                    my ($pre, $post) = ();
                    if ($ai_var->{sample}->{$sample}->{PLOIDY} eq '1') {
                        $pre = "hap";
                    }
                    else {
                        $pre = "hom";
                    }

                    $ai_var->{sample}->{$sample}->{zygosity} =
                      $pre . '-' . $cur_post;
                }
                elsif (!exists $locus_strand_index->{'3'}) {
                    my $pair_ai = ($ai eq '1') ? '2' : '1';
                    if (!exists $locus_strand_index->{$pair_ai}) {
                        if ( $ai_var->{sample}->{$sample}->{PLOIDY} eq '1') {
                            $ai_var->{sample}->{$sample}->{zygosity} = 'hap-'.$cur_post;
                        }
                        else {
                            # not complete locus?
                            # set default to '.'
                            carp "[Critical Warning] may not complete locus $rlg->[0]->{locus}, set zygosity to '.' as unknown.\n";
                            $ai_var->{sample}->{$sample}->{zygosity} = '.';
                        }
                    }
                    else {
                        # select pair var
                        my $pair_var;
                        foreach my $pair_var_idx (
                            @{ $locus_strand_index->{$pair_ai} } )
                        {
                            # select the first overlapped as pair_var
                            if ( $ai_var->{begin} <
                                    $rlg->[$pair_var_idx]->{end}
                                and $ai_var->{end} >
                                $rlg->[$pair_var_idx]->{begin} )
                            {
                                $pair_var = $rlg->[$pair_var_idx];
                                last;
                            }
                        }

                        if ( !defined $pair_var ) {

                            # not found overlapped pair var
                            # assume overlapped with reference
                            $ai_var->{sample}->{$sample}->{zygosity} =
                              ( $cur_post eq 'ref' ) ? 'hom-ref' : 'het-ref';
                        }
                        else {

                            if ( $pair_var->{variantSequence} eq '?' ) {
                                $ai_var->{sample}->{$sample}->{zygosity} =
                                  'half-' . $cur_post;
                            }
                            else {
                                if ( $pair_var->{callerVarType} =~ /ref/ ) {
                                    $ai_var->{sample}->{$sample}->{zygosity} =
                                      ( $cur_post eq 'ref' )
                                      ? 'hom-ref'
                                      : 'het-ref';
                                }
                                else {
                                    $ai_var->{sample}->{$sample}->{zygosity} =
                                      'het-alt';
                                }
                            }
                        }
                    }
                }
                else {
                    # multiple ploidy will always with unknown zygosity
                    # except for ai_all and no-call
                    $ai_var->{sample}->{$sample}->{zygosity} = '.';
                }
            }
        }
    }

    for (my $i = $#$rlg; $i >= 0; $i--) {
        # delete long ref call for current protocol
        if ($rlg->[$i]->{referenceSequence} eq '='
                and $rlg->[$i]->{variantSequence} eq '=') {
            splice( @$rlg, $i, 1);
        }
    }

    return $rlg;
}

# combind hom-alt mutation
sub check_combin_group {
    my ( $rlg, $rlsi ) = @_;
    return ($rlg, $rlsi) if ( 2 != scalar keys %$rlsi or 2 != scalar @$rlg );
    return ($rlg, $rlsi) if ( !exists $rlsi->{1} or !exists $rlsi->{2} );
    my $v1 = $rlg->[0];
    my $v2 = $rlg->[1];
    my @samples = keys %{$v1->{sample}};
    if (    $v1->{begin} eq $v2->{begin}
        and $v1->{end} eq $v2->{end}
        and $v1->{variantSequence} eq $v2->{variantSequence} )
    {
        foreach my $sam (@samples) {
            if ( $v1->{sample}->{$sam}->{filterTag} ne $v2->{sample}->{$sam}->{filterTag} ) {
                if (   ( $v1->{sample}->{$sam}->{filterTag} ne 'PASS' and $v1->{sample}->{$sam}->{filterTag} ne '.' )
                    or ( $v2->{sample}->{$sam}->{filterTag} ne 'PASS' and $v2->{sample}->{$sam}->{filterTag} ne '.' ) )
                {
                    $rlg->[0]->{sample}->{$sam}->{filterTag} = 'DUBIOUS';
                }
            }

            if ( $v1->{sample}->{$sam}->{VF} ne $v2->{sample}->{$sam}->{VF} ) {
                $rlg->[0]->{sample}->{$sam}->{VF} .= "|".$rlg->[1]->{sample}->{$sam}->{VF};
            }

            $rlg->[0]->{sample}->{$sam}->{AI} = 0;
        }
        delete $rlg->[1];
        $rlsi = { 0 => [0] };
    }

    return ($rlg, $rlsi);
}

sub tsv_var_parser {
    my ($ritms, $sample) = @_;

    my $var = {
        locus             => $$ritms[0],
        chr               => $$ritms[3],
        begin             => $$ritms[4],
        end               => $$ritms[5],
        callerVarType     => $$ritms[6],
        referenceSequence => $$ritms[7],
        variantSequence   => $$ritms[8],
    };

    my ($AI, $PB, $VF, $PLOIDY);
    $PB = ($$ritms[12] eq '') ? '.' : $$ritms[12];
    $VF = ($$ritms[11] eq '') ? '.' : $$ritms[11];
    $AI = ($$ritms[2] eq '')  ? '.' : $$ritms[2];
    $AI = ($AI eq 'all') ? 0 : $AI;
    $PLOIDY = $$ritms[1];

    $var->{sample}->{$sample} = {
        PLOIDY   => $PLOIDY,
        VF       => $VF,
        AI       => $AI,
        PB       => $PB,
    };
    $var->{sample}->{$sample} = give_filter( $var->{sample}->{$sample} );

    return $var;
}

# Current InExcel Strategy 
# will depend on autoInterpRules
# defined in config
# Or default rules will be implemented

sub check_inExcel {
    my $rec_ent = shift;
    if (
        exists $kickoutExcel_function{ $rec_ent->{Function} }
        and (  $rec_ent->{autoInterp} =~ /Benign/
            or $rec_ent->{autoInterp} eq "Unknown" )
      )
    {
        if ( $rec_ent->{Function} eq "intron" ) {
            return 1 if ($config{intron_edge_remain} <= 0);
            while ($rec_ent->{cHGVS} =~ /\d+[\+\-](\d+)/g) {
                if ($1 < $config{intron_edge_remain}) {
                    return 1;
                }
            }
        }
        return 0;
    }

    return 1;
}

# AutoInterp status comes from the following options:
# 1. whether related to a DM mutation record
# 2. whether a high quality calling
# 3. whether with a deleterious function
# 4. whether have lower allele frequency than threshold
# if all are yes: then "Disease Causing Mutation",
# if no for all without consider quality, then "Benign",
# if 1,4 or 3,4 are yes: then "Possibly Deleterious"
# if only 2,1(and/or)3, only 2,4 are yes, then "Possibly Benign",
# other case, then "Unknown"

sub AutoInterp {
    my $rec_ent = shift;

    my $Phyno_pred_opt = 0;
    my $Func_pred_opt = 0;
    my $freq_opt   = 1;
    my $filter_opt = 0;

    # phyno_pred hgmd and GaP
    if (
           $rec_ent->{hgmdPred} =~ /DM/
        or $rec_ent->{hgmdPred} =~ /DFP/
        or $rec_ent->{ClinSignificant} =~ /pathogenic/i )
    {
        $Phyno_pred_opt = 1;
    }

    # filter
    if ($rec_ent->{Filter} eq 'PASS' or $rec_ent->{Filter} eq '.') {
        $filter_opt = 1;
    }

    # func
    if (
           exists $deleterious_func{ $rec_ent->{Function} }
        or exists $possible_deleterious_func{ $rec_ent->{Function} }
        or $rec_ent->{ensCondelPred} eq 'deleterious'
      )
    {
        $Func_pred_opt = 1;
    }

    # freq
    foreach my $freq_k ( keys %autoInterp_freq_threshold ) {
        if (    exists $rec_ent->{$freq_k}
            and $rec_ent->{$freq_k} ne '.'
            and $rec_ent->{$freq_k} ne 'error'
            and $rec_ent->{$freq_k} > $autoInterp_freq_threshold{$freq_k} )
        {
            $freq_opt = 0;
            last;
        }
    }

    # if all are yes: then "Certain" (Disease Mutation),
    if (    $Phyno_pred_opt == 1
        and $filter_opt == 1
        and $Func_pred_opt == 1
        and $freq_opt == 1 )
    {
        return "Certain";
    }

    # if no for all without consider quality, then "Benign",
    elsif ( $Phyno_pred_opt == 0 and $Func_pred_opt == 0 and $freq_opt == 0 ) {
        return "Benign";
    }

    # if 1,4 are yes: then "Likely Deleterious"
    elsif ( $Phyno_pred_opt == 1 and $freq_opt == 1 )
    {
        return "Likely Deleterious";
    }

    # if 3,4 are yes: then "VUS"
    elsif ( $Func_pred_opt == 1 and $freq_opt == 1 )
    {
        return "VUS";
    }

    # if only 2,1(and/or)3, only 2,4 are yes, then "Likely Benign",
    elsif ( $filter_opt == 1
        and
        ( ( $Phyno_pred_opt == 1 or $Func_pred_opt == 1 ) xor $freq_opt == 1 ) )
    {
        return "Likely Benign";
    }

    # other case, then "Unknown"
    else {
        return "Unknown";
    }
}

sub checkXPAR {
    my ($begin, $end) = @_;
    foreach my $rxpar (@xpar_reg) {
        if ($begin < $rxpar->[1] and $end > $rxpar->[0]) {
            return 0;
        }
    }
    return 1;
}


sub checkYPAR {
    my ($begin, $end) = @_;
    foreach my $rypar (@ypar_reg) {
        if ($begin < $rypar->[1] and $end > $rypar->[0]) {
            return 0;
        }
    }
    return 1;
}

sub uniform_chrVarOut {
    my ($varanno, $sampleid) = @_;
    my $varout = {};
    my $chrvar = $varanno->{var};
    my $sample_h = $chrvar->{sample}->{$sampleid};
    my $ExAC_pop   = $config{ExAC_POP};
    my $TGP_pop = $config{TGP_POP};
    my $GAD_pop = 'EAS'; # api only support EAS

    # must have keys in var
    @$varout{qw(
        Chr Start Stop Ref 
        VarType Call
    )} = @$chrvar{qw(
        chr pos end ref
        guess alt
    )};

    my @tmpids = (
    qw(
        Flank MapLoc RepeatTag
        dbsnpNote dbsnpAC dbsnpAF dbsnpMTag
        TGP_pop_AF TGP_AF ESP6500AC ESP6500AF
        ExAC_Filter ExAC_pop_HASC ExAC_pop_AC
        ExAC_pop_AF ExAC_HASC ExAC_AC ExAC_AF
        GAD_pop_HASC GAD_pop_AC
        GAD_pop_AF GAD_HASC GAD_AC GAD_AF
        rsID PVFD_RULE
        PVFDhomaC PVFDhomaF PVFD_AC PVFD_AF
        PanelID PanelAC PanelAF 
        phyloPpr phyloPve phyloPpm
        cosmicName cosmicSite cosmicHis cosmicStat cosmicID
        ClinSSR ClinACC ClinRevStat ClinSignificant
    ));
    @$varout{@tmpids} = ('.') x (scalar @tmpids);
    
    if ($chrvar->{guess} ne "annotation-error") {
        # must have if annotation OK
        @$varout{qw(
            Flank MapLoc
        )} = @$chrvar{qw(
            flanks cytoBand
        )};
    }

    # for single point only
    if (exists $chrvar->{phyloPpr} and defined $chrvar->{phyloPpr}) {
        @$varout{qw(
            phyloPpr phyloPve phyloPpm
        )} = @$chrvar{qw(
            phyloPpr phyloPve phyloPpm
        )};
    }

    # keys for sample 
    $varout->{SampleID} = $sampleid;
    if ($sample_list) {
        $varout->{FamID} =
          ( exists $sample_info->{$sampleid} )
          ? $sample_info->{$sampleid}->{fam}
          : '.';
        $varout->{Sex} =
          ( exists $sample_info->{$sampleid} )
          ? $sample_info->{$sampleid}->{sex}
          : '.';
    }
    
    @$varout{qw(
        Zygosity PhasedGID A.Index Filter
    )} = @$sample_h{qw(
        zygosity PB AI filterTag
    )};

    # fix zygosity info for male's sex chromosome
    if (
            exists $varout->{Sex}
        and $varout->{Sex} eq 'M'
        and (
            (
                $varout->{Chr} =~ /^Y/i
                and checkYPAR( $varout->{Start}, $varout->{Stop} )
            )
            or ( $varout->{Chr} =~ /^X/i
                and checkXPAR( $varout->{Start}, $varout->{Stop} ) )
        )
      )
    {
        $varout->{Zygosity} =~ s/^hom-/hem-/;
    }

    $varout->{NbGID} = ((exists $sample_h->{NB}) ? $sample_h->{NB} : '.');
    $varout->{"A.Depth"} = ((exists $sample_h->{AD}) ? $sample_h->{AD} : '.');
    $varout->{"A.Ratio"} = ((exists $sample_h->{AR}) ? $sample_h->{AR} : '.');

    if (exists $chrvar->{reptag}) {
        $varout->{RepeatTag} = $chrvar->{reptag};
    }

    # keys for AF info
    if ( exists $chrvar->{dbsnp}
        and 0 < scalar keys %{ $chrvar->{dbsnp} } )
    {
        my @rs_ids = sort keys %{ $chrvar->{dbsnp} };
        my ( @dbsnpNotes, @dbsnpACs, @dbsnpAFs, @dbsnpMTags ) = ();
        foreach my $rs (@rs_ids) {
            my $cur_rs_h = $chrvar->{dbsnp}->{$rs};
            push(
                @dbsnpNotes,
                (
                    ( $cur_rs_h->{exception} eq "." ) ? $cur_rs_h->{bitfields}
                    : (
                        ( $cur_rs_h->{bitfields} eq "." )
                        ? $cur_rs_h->{exception}
                        : $cur_rs_h->{exception} . ";" . $cur_rs_h->{bitfields}
                    )
                )
            );
            push( @dbsnpACs, $cur_rs_h->{AN} );
            push( @dbsnpAFs, $cur_rs_h->{AF} );
            push(
                @dbsnpMTags,
                (
                    ( $cur_rs_h->{weight} eq "1" )
                    ? "U"
                    : "M"
                )
            );
        }
        $varout->{rsID}      = subsame(@rs_ids);
        $varout->{dbsnpNote} = subsame(@dbsnpNotes);
        $varout->{dbsnpAC}   = subsame(@dbsnpACs);
        $varout->{dbsnpAF}   = subsame(@dbsnpAFs);
        $varout->{dbsnpMTag} = subsame(@dbsnpMTags);
    }

    if (    exists $chrvar->{tgp}
        and exists $chrvar->{tgp}->{AF} )
    {
        $varout->{TGP_AF}     = $chrvar->{tgp}->{AF};
        $varout->{TGP_pop_AF} = $chrvar->{tgp}->{ "AF_" . $TGP_pop }
          if ( exists $chrvar->{tgp}->{ "AF_" . $TGP_pop } );
        $varout->{TGP_HSAC} = $chrvar->{tgp}->{HSN}
          if ( exists $chrvar->{tgp}->{HSN} );
    }

    if (    exists $chrvar->{esp6500}
        and exists $chrvar->{esp6500}->{AF} )
    {
        $varout->{ESP6500AC} = $chrvar->{esp6500}->{AN};
        $varout->{ESP6500AF} = $chrvar->{esp6500}->{AF};
    }

    if (    exists $chrvar->{exac}
        and exists $chrvar->{exac}->{AF} )
    {
        my $exac_hash = $chrvar->{exac};

        $varout->{ExAC_Filter} = $exac_hash->{Filter};
        $varout->{ExAC_pop_HASC} = $exac_hash->{ "HASC_" . $ExAC_pop }
          if ( exists $exac_hash->{ "HASC_" . $ExAC_pop } );
        $varout->{ExAC_pop_AC} = $exac_hash->{ "AN_" . $ExAC_pop }
          if ( exists $exac_hash->{ "AN_" . $ExAC_pop } );
        $varout->{ExAC_pop_AF} = $exac_hash->{ "AF_" . $ExAC_pop }
          if ( exists $exac_hash->{ "AF_" . $ExAC_pop } );
        $varout->{ExAC_HASC} = $exac_hash->{HASC}
          if ( exists $exac_hash->{HASC} );
        $varout->{ExAC_AC} = $exac_hash->{AN}
          if ( exists $exac_hash->{AN} );
        $varout->{ExAC_AF} = $exac_hash->{AF}
          if ( exists $exac_hash->{AF} );
    }

    if (    exists $chrvar->{gnomAD}
        and exists $chrvar->{gnomAD}->{AF} )
    {
        my $gnomAD_hash = $chrvar->{gnomAD};

        $varout->{GAD_pop_HASC} = $gnomAD_hash->{ "HASC_" . $GAD_pop }
          if ( exists $gnomAD_hash->{ "HASC_" . $GAD_pop } );
        $varout->{GAD_pop_AC} = $gnomAD_hash->{ "AN_" . $GAD_pop }
          if ( exists $gnomAD_hash->{ "AN_" . $GAD_pop } );
        $varout->{GAD_pop_AF} = $gnomAD_hash->{ "AF_" . $GAD_pop }
          if ( exists $gnomAD_hash->{ "AF_" . $GAD_pop } );
        $varout->{GAD_HASC} = $gnomAD_hash->{HASC}
          if ( exists $gnomAD_hash->{HASC} );
        $varout->{GAD_AC} = $gnomAD_hash->{AN}
          if ( exists $gnomAD_hash->{AN} );
        $varout->{GAD_AF} = $gnomAD_hash->{AF}
          if ( exists $gnomAD_hash->{AF} );
    }

    if (exists $chrvar->{PanelSQL_Fail}) {
        $varout->{PanelAF} = 'error';
    }

    if (    exists $chrvar->{panelctrl}
        and exists $chrvar->{panelctrl}->{AF} )
    {
        $varout->{PanelID} = $config{Panel_ID};
        $varout->{PanelAC} = $chrvar->{panelctrl}->{AN};
        $varout->{PanelAF} = $chrvar->{panelctrl}->{AF};
    }

    # clinVar
    if (exists $chrvar->{clinVar_Fail}) {
        $varout->{ClinSignificant} = 'error';
    }
    if ( exists $chrvar->{clinVar} and exists $chrvar->{clinVar}->{SSR} ) {
        @$varout{ qw(ClinSSR ClinACC ClinRevStat ClinSignificant) } =
          @{ $chrvar->{clinVar} }{ qw(SSR CLNACC CLNREVSTAT CLNSIG) };
    }

    # cosmic
    if (    exists $chrvar->{cosmic}
        and exists $chrvar->{cosmic}->{mutID} )
    {
        @$varout{
            qw(
              cosmicName cosmicSite cosmicHis cosmicStat cosmicID
              )
          }
          = @{ $chrvar->{cosmic} }{
            qw(
              mutName site histology status mutID
              )
          };
    }

    return $varout;
}

sub uniform_trVarOut {
    my ($varanno, $varout, $tr) = @_;
    
    my $new_var_out = {%$varout};
    my @ids = (
    qw(
        mimGID mimStat mimPIDs mimInhs
        cgdCond cgdInhs cgdACond cgdManCat cgdIntCat cgdRef
        TestCode EntrezGeneID geneSym
        FuncRegion ExIn_ID pfamName pfamId Function
        Transcript Protein Strand Primary cHGVS pHGVS
        CodonChange PolarChange MutationName
        ensSIFT ensPP2hvar ensPP2hdiv ensCondel ensCondelPred
        hgmdID hgmdQMode hgmdName hgmdDis hgmdPMID hgmdPred
    ));
    @$new_var_out{@ids} = ('.') x (scalar @ids);

    if ($tr ne "") {

        my $curtr = $varanno->{trInfo}->{$tr};
        $new_var_out->{Transcript} = $tr;

        # must defined keys
        @$new_var_out{qw(
            EntrezGeneID geneSym TestCode
            FuncRegion ExIn_ID Function
            Protein Strand Primary
            cHGVS MutationName
        )} = @$curtr{qw(
            geneId geneSym TCD
            r exin func 
            prot strd major
            c trVarName
        )};


        my $geneVarNameRef;
        my $custom_VarName;
        if (exists $config{VarNameList}) {
            if (exists $VarNameList{$tr}) {
                $geneVarNameRef = $VarNameList{$tr};
            }
            elsif (exists $VarNameList{$$curtr{geneSym}}) {
                $geneVarNameRef = $VarNameList{$$curtr{geneSym}};
            }
            elsif (exists $VarNameList{$$curtr{geneId}}) {
                $geneVarNameRef = $VarNameList{$$curtr{geneId}};
            }

            if (defined $geneVarNameRef) {
                if (exists $geneVarNameRef->{$$curtr{c}}) {
                    $custom_VarName = $geneVarNameRef->{$$curtr{c}};
                }
                else {
                    my $tmp = $$curtr{c};

                    if ($$curtr{c} =~ /\-u/) {
                        $tmp =~ s/\-(\d+)\-u(\d+)/"-".($1+$2)/eg;
                        $tmp =~ s/\-u/\-/g; # for those no 5'utr genes
                    }

                    if ($$curtr{c} =~ /\+d/) {
                        $tmp =~ s/\*(\d+)\+d(\d+)/"*".($1+$2)/eg;
                        $tmp =~ s/\+d/\*/g; # for those no 3'utr genes
                    }

                    if ($$curtr{c} =~ /dup/) {
                        $tmp =~ s/dup.+$/dup/; # for other format of dup
                    }

                    if ($$curtr{c} =~ /del/) {
                        $tmp =~ s/del.+$/del/;
                    }

                    if (defined $geneVarNameRef and exists $geneVarNameRef->{$tmp}) {
                        $custom_VarName = $geneVarNameRef->{$tmp};
                    }
                }
            }
        }

        if (exists $curtr->{standard_cHGVS} or exists $curtr->{alt_cHGVS}) {
            $new_var_out->{cHGVS} .= ' (';
            $new_var_out->{cHGVS} .= 'std: '.$curtr->{standard_cHGVS}.' ' if (exists $curtr->{standard_cHGVS});
            $new_var_out->{cHGVS} .= 'alt: '.$curtr->{alt_cHGVS}.' ' if (exists $curtr->{alt_cHGVS});
            $new_var_out->{cHGVS} .= ')';

            if (defined $geneVarNameRef and !defined $custom_VarName) {
                if (exists $geneVarNameRef->{$curtr->{standard_cHGVS}}) {
                    $custom_VarName = $geneVarNameRef->{$curtr->{standard_cHGVS}};
                }
                elsif (exists $geneVarNameRef->{$curtr->{alt_cHGVS}}) {
                    $custom_VarName = $geneVarNameRef->{$curtr->{alt_cHGVS}};
                }
            }
        }

        if (exists $curtr->{p}) {
            $new_var_out->{pHGVS} = $curtr->{p};
            if (    defined $geneVarNameRef
                and !defined $custom_VarName
                and exists $geneVarNameRef->{ $curtr->{p} } )
            {
                $custom_VarName = $geneVarNameRef->{ $curtr->{p} };
            }
        }

        if (exists $curtr->{standard_pHGVS} or exists $curtr->{alt_pHGVS}) {
            $new_var_out->{pHGVS} .= ' (';
            $new_var_out->{pHGVS} .= 'std: '.$curtr->{standard_pHGVS}.' ' if (exists $curtr->{standard_pHGVS});
            $new_var_out->{pHGVS} .= 'alt: '.$curtr->{alt_pHGVS}.' ' if (exists $curtr->{alt_pHGVS});
            $new_var_out->{pHGVS} .= ')';

            if (defined $geneVarNameRef and !defined $custom_VarName) {
                if (exists $curtr->{strandard_pHGVS} and exists $geneVarNameRef->{ $curtr->{standard_pHGVS} }) {
                    $custom_VarName = $geneVarNameRef->{ $curtr->{standard_pHGVS} };
                }
                elsif (exists $curtr->{alt_pHGVS} and exists $geneVarNameRef->{ $curtr->{alt_pHGVS} } ) {
                    $custom_VarName = $geneVarNameRef->{ $curtr->{alt_pHGVS} };
                }
            }
        }

        if (exists $curtr->{p3} and $curtr->{p3} ne $curtr->{p}) {
            $new_var_out->{pHGVS} .= ' | '.$curtr->{p3};
            if (defined $geneVarNameRef and !defined $custom_VarName and exists $geneVarNameRef->{$$curtr{p3}}) {
                $custom_VarName = $geneVarNameRef->{$$curtr{p3}};
            }
            if (exists $curtr->{standard_p3} or exists $curtr->{alt_p3}) {
                $new_var_out->{pHGVS} .= ' (';
                $new_var_out->{pHGVS} .= 'std: '.$curtr->{standard_p3}.' ' if (exists $curtr->{standard_p3});
                $new_var_out->{pHGVS} .= 'alt: '.$curtr->{alt_p3}.' ' if (exists $curtr->{alt_p3});
                $new_var_out->{pHGVS} .= ')';
                if (defined $geneVarNameRef and !defined $custom_VarName) {
                    if (exists $geneVarNameRef->{ $curtr->{standard_p3} }) {
                        $custom_VarName = $geneVarNameRef->{ $curtr->{standard_p3} };
                    }
                    elsif (exists $geneVarNameRef->{ $curtr->{alt_p3} } ) {
                        $custom_VarName = $geneVarNameRef->{ $curtr->{alt_p3} };
                    }
                }
            }
        }

        if (defined $custom_VarName) {
            $new_var_out->{MutationName} .= ' / '. $custom_VarName;
        }
        $new_var_out->{CodonChange} = $curtr->{cc} if ( exists $curtr->{cc} );
        $new_var_out->{PolarChange} = $curtr->{polar}
          if ( exists $curtr->{polar} );

        # pfam
        if ( exists $curtr->{pfamId} and defined $curtr->{pfamId} ) {
            @$new_var_out{qw(
                pfamId pfamName 
            )} = @$curtr{qw(
                pfamId pfamName
            )};
        }

        # sift condel
        if ( exists $curtr->{condelPred} and defined $curtr->{condelPred}) {
            @$new_var_out{qw(
                ensSIFT ensPP2hvar ensPP2hdiv ensCondel ensCondelPred
            )} = @$curtr{qw(
                siftScore pp2varScore pp2divScore 
                condelScore condelPred
            )};
        }


        if (exists $curtr->{GeneSQL_Fail}) {
            $new_var_out->{mimInhs} = 'error';
            $new_var_out->{cgdRef}  = 'error';
        }

        # mim
        if (exists $curtr->{omim} and 0 < scalar keys %{$curtr->{omim}}) {
            my @genemims = sort keys %{$curtr->{omim}};
            my (@genestat, @disomims, @inhs) = ();
            foreach my $gm (@genemims) {
                $curtr->{omim}->{$gm}->{genestat} =~ tr/|/;/;
                $curtr->{omim}->{$gm}->{disomims} =~ tr/|/;/;
                $curtr->{omim}->{$gm}->{inhs} =~ tr/|/;/;

                push (@genestat, $curtr->{omim}->{$gm}->{genestat});
                push (@disomims, $curtr->{omim}->{$gm}->{disomims});
                push (@inhs,     $curtr->{omim}->{$gm}->{inhs});
            }

            @$new_var_out{qw(
                mimGID mimStat mimPIDs mimInhs
            )} = (
                join('|', @genemims),
                subsame(@genestat),
                subsame(@disomims),
                subsame(@inhs)
            );
        }

        # cgd
        if (exists $curtr->{cgd} and exists $curtr->{cgd}->{geneid}) {
            $new_var_out->{cgdRef} = trim_pmID($curtr->{cgd}->{references});
            @$new_var_out{qw(
                cgdCond cgdInhs cgdACond 
                cgdManCat cgdIntCat
            )} = @{$curtr->{cgd}}{qw(
                cond inhs allelic_cond
                manifest_cat intervent_cat
            )};
        }

        if (exists $curtr->{LocalHGMD_Fail}) {
            $new_var_out->{hgmdPred} = 'error';
        }

        # hgmd
        if (exists $curtr->{hgmd} and 0 < (scalar @{$curtr->{hgmd}})) {
            my ( @hgmdID, @hgmdQMode, @hgmdName,
                @hgmdDis, @hgmdPMID, @hgmdPred ) = ();
            my @all_ret = sort { $a->{id} cmp $b->{id} } @{ $curtr->{hgmd} };
            foreach my $ret1 (@all_ret) {
                push( @hgmdID,    $ret1->{id} );
                push( @hgmdQMode, $ret1->{querylevel} )
                  if ( exists $ret1->{querylevel} );
                push( @hgmdName, $ret1->{mutation_name} );
                push( @hgmdDis,  $ret1->{disease} );
                push( @hgmdPMID, $ret1->{pmid} );
                push( @hgmdPred, $ret1->{pred} );
            }
            $new_var_out->{hgmdID}    = subsame(@hgmdID);
            $new_var_out->{hgmdQMode} = subsame(@hgmdQMode)
              if ( defined $hgmdQMode[0] );
            $new_var_out->{hgmdName} = subsame(@hgmdName);
            $new_var_out->{hgmdDis}  = subsame(@hgmdDis);
            $new_var_out->{hgmdPMID} = trim_pmID(subsame(@hgmdPMID));
            $new_var_out->{hgmdPred} = subsame(@hgmdPred);
        }

    }
    else {
        $new_var_out->{MutationName} = $varanno->{var}->{varName}
          if ( exists $varanno->{var}->{varName} );

    }

    return $new_var_out;
}

sub trim_pmID {
    my $pmID = shift;
    my %has_pmID = ();
    while ($pmID =~ /(\d+)/g) {
        $has_pmID{$1} = 1;
    }
    my $rejoin = join('|', sort {$b<=>$a} keys %has_pmID);
    return $rejoin;
}

sub subsame {
    my ($val_ori, @string_array) = @_;
    my $same_opt = 1;
    foreach (@string_array) {
        if ($val_ori ne $_) {
            $same_opt = 0;
            last;
        }
    }
    if (0 == scalar(@string_array) or $same_opt == 1) {
        return $val_ori;
    }
    else {
        return join('|', $val_ori, @string_array);
    }
}

# count while output
sub print_cache {
    # already same chromosome, only need to sort position
    my @sort_anno = sort {
             $a->{var}->{pos} <=> $b->{var}->{pos}
          or $a->{var}->{end} <=> $b->{var}->{end}
    } @anno_cache;
    
    foreach my $var_ent (@sort_anno) {
        
        foreach my $sample (sort keys %{$var_ent->{var}->{sample}}) {
            my $curSample_h = $var_ent->{var}->{sample}->{$sample};
            my $invar_allele = 1;
            if ($curSample_h->{AI} == 0 and $curSample_h->{PLOIDY} == 2) {
                $invar_allele = 2;
            }

            if ( $var_ent->{var}->{guess} eq 'snv' ) {
                if (
                    judge_titv(
                        $var_ent->{var}->{ref}, $var_ent->{var}->{alt}
                    )
                  )
                {
                    $ti{$sample} += $invar_allele;
                }
                else {
                    $tv{$sample} += $invar_allele;
                }
            }

            if ($var_ent->{var}->{chr} =~ /^[xy]/i) {
                # sex chromosome
                if ( $var_ent->{var}->{guess} ne 'no-call' ) {
                    $sex_varcount{$sample}{total} += $invar_allele;

                    $sex_varcount{$sample}{het} += $invar_allele
                      if ( $curSample_h->{zygosity} =~ /het/ );
                }
            }

            if ( $var_ent->{var}->{guess} eq 'no-call' ) {
                $no_call_var_count{$sample} += $invar_allele;
            }
            else {
                $var_count{$sample}{total} += $invar_allele;
                $var_count{$sample}{varType}{ $var_ent->{var}->{guess} } +=
                  $invar_allele;
                $var_count{$sample}{zygosity}{ $curSample_h->{zygosity} } +=
                  $invar_allele;
                $var_count{$sample}{filter}{ $curSample_h->{filterTag} } +=
                  $invar_allele;
            }

            my $arbitrary_remain;
            if (defined $coreFuncRegion_h) {
                if (
                    $coreFuncRegion_h->checkIn(
                        $var_ent->{var}->{chr}, $var_ent->{var}->{pos},
                        $var_ent->{var}->{end}
                    )
                )
                {
                    $arbitrary_remain = 1;
                }
            }

            
            my $cur_var_out = uniform_chrVarOut($var_ent, $sample);
            if ( exists $var_ent->{trInfo} ) {
                foreach my $tr (sort keys %{$var_ent->{trInfo}}) {
                    my $tr_var_out =
                      uniform_trVarOut( $var_ent, $cur_var_out, $tr );
                    $tr_var_out->{autoInterp} = AutoInterp($tr_var_out);

                    # must have autoInterp first
                    if (defined $arbitrary_remain and $arbitrary_remain == 1) {
                        $tr_var_out->{InExcel} = 1;
                    }
                    else {
                        $tr_var_out->{InExcel} = check_inExcel($tr_var_out);
                    }

                    # only count the major transcript for called-allele
                    # to avoid multiple transcript dup
                    if (    $var_ent->{var}->{guess} ne 'no-call'
                        and $var_ent->{var}->{varName} eq
                        $tr_var_out->{MutationName} )
                    {
                        $annotation_count{$sample}{autoInterp}
                          { $tr_var_out->{autoInterp} } += $invar_allele;
                        $annotation_count{$sample}{function}
                          { $tr_var_out->{Function} } += $invar_allele;
                        $annotation_count{$sample}{InExcel}
                          { $tr_var_out->{InExcel} } += $invar_allele;
                        $annotation_count{$sample}{total} += $invar_allele;
                    }

                    print STDERR "uniformed var: ".Dumper($tr_var_out) if (defined $debug);

                    print $OutFp join("\t", normalise_null(@$tr_var_out{@out_header_keys}))."\n";
                }
            }
            else {
                $cur_var_out = uniform_trVarOut( $var_ent, $cur_var_out, "" );
                $cur_var_out->{autoInterp} = AutoInterp($cur_var_out);

                # must have autoInterp first
                if ( $var_ent->{var}->{guess} ne 'no-call' ) {
                    if (defined $arbitrary_remain and $arbitrary_remain == 1) {
                        $cur_var_out->{InExcel} = 1;
                    }
                    else {
                        $cur_var_out->{InExcel} = check_inExcel($cur_var_out);
                    }
                    $annotation_count{$sample}{autoInterp}
                      { $cur_var_out->{autoInterp} } += $invar_allele;
                    $annotation_count{$sample}{function}
                      { $cur_var_out->{Function} } += $invar_allele;
                    $annotation_count{$sample}{InExcel}
                      { $cur_var_out->{InExcel} } += $invar_allele;
                    $annotation_count{$sample}{total} += $invar_allele;
                }
                print STDERR "uniformed var: ".Dumper($cur_var_out) if (defined $debug);
                print $OutFp join("\t", normalise_null(@$cur_var_out{@out_header_keys}))."\n";
            }
        }
    }
}

# print the header line to the Output
sub print_header {

    @out_header_keys = (); # clean output header keys
    my $TGP_pop    = $config{TGP_POP};
    my $ExAC_pop   = $config{ExAC_POP};
    my $GAD_pop = 'EAS'; #current api only support EAS population

    # =========== Location group ========= #
    my @chrpos_header = ( "Chr", "Start", "Stop", "MapLoc" );
    my $inExcel = "InExcel";    # firm specified.
    my @sample_header;
    if ($sample_list) {
        @sample_header = ( "FamID", "SampleID", "Sex" );
    }
    else {
        @sample_header = ("SampleID");
    }

    # ========== Genotyping Group ========= #
    my @geno_header = ();
    push( @geno_header,
        "NbGID",     "Ref",      "VarType",   "Call",
        "Flank",     "Zygosity", "A.Depth",   "A.Ratio",
        "PhasedGID", "A.Index",  "RepeatTag", "Filter" );

    push( @out_header_keys,
        @chrpos_header, $inExcel, @sample_header, @geno_header );

    # ========== Gene Information ========= #
    my @gene_header = (
        "MIM Gene ID",
        "MIM Stat",
        "MIM Pheno IDs",
        "MIM Inheritance",
        "CGD Condition",
        "CGD Inheritance",
        "CGD Allelic Condition",
        "CGD Manifestation Categories",
        "CGD Intervention Categories",
        "CGD References",
        "TestCode",
        "EntrezGeneID",
        "Gene Symbol"
    );

    push( @out_header_keys, qw{
            mimGID mimStat mimPIDs mimInhs
            cgdCond cgdInhs cgdACond cgdManCat cgdIntCat cgdRef
            TestCode EntrezGeneID geneSym
        } );

    # ======== Function Information ======== #
    my @func_header =
      ( "FuncRegion", "ExIn_ID", "pfamName", "pfamId", "Function" );

    # ========== HGVS Information ========== #
    my @hgvs_header = (
        "Transcript", "Protein", "Strand",      "Primary",
        "cHGVS",      "pHGVS",   "CodonChange", "PolarChange",
        "MutationName"
    );

    push( @out_header_keys, @func_header, @hgvs_header );

    # ========= Public Frequency DB ======== #
    my @freq_header = ();

    # dbSNP
    push( @freq_header,
        "dbSNP Note",
        "dbSNP Allele Count",
        "dbSNP Allele Freq",
        "dbSNP Mul Tag" );

    # 1000 genomes
    push( @freq_header, "1000G $TGP_pop AF", "1000G AF" );

    # ESP6500
    push( @freq_header, "ESP6500 AC", "ESP6500 AF" );

    # ExAC
    push( @freq_header, 
        "ExAC Filter", 
        "ExAC $ExAC_pop HomoAlt Count",
        "ExAC $ExAC_pop AC",
        "ExAC $ExAC_pop AF",
        "ExAC HomoAlt Count",
        "ExAC AC",
        "ExAC AF" );

    # GAD
    push( @freq_header, 
        "GAD $GAD_pop HomoAlt Count",
        "GAD $GAD_pop AC",
        "GAD $GAD_pop AF",
        "GAD HomoAlt Count",
        "GAD AC",
        "GAD AF" );

    # rsID from dbSNP
    push( @freq_header, "rsID" );

    push( @out_header_keys, qw{
            dbsnpNote dbsnpAC dbsnpAF dbsnpMTag
            TGP_pop_AF TGP_AF ESP6500AC ESP6500AF
            ExAC_Filter ExAC_pop_HASC ExAC_pop_AC
            ExAC_pop_AF ExAC_HASC ExAC_AC ExAC_AF
            GAD_pop_HASC GAD_pop_AC
            GAD_pop_AF GAD_HASC GAD_AC GAD_AF
            rsID
        } );

    # ========= Panel Frequency DB ========= #
    my @panel_header = ( "Panel ID", "Panel AlleleCount", "Panel AlleleFreq" );

    push(
        @out_header_keys, qw{
          PanelID PanelAC PanelAF
          }
    );

    # ====== Mutation Effect Prediction ==== #

    # PhyloP scores
    my @phyloP_header =
      ( "PhyloP Primates", "PhyloP Vertebrates", "PhyloP Placental Mammals" );

    push( @out_header_keys, qw{
            phyloPpr phyloPve phyloPpm
        } );

    # Condel prediction
    my @condel_header = (
        "Ens SIFT Score",
        "Ens Polyphen2HumVar Score",
        "Ens Polyphen2HumDiv Score",
        "Ens Condel Score",
        "Ens Condel Pred"
    );

    # Cosmic prediction
    my @cosmic_header = (
        "Cosmic MutName",
        "Cosmic Site",
        "Cosmic Histology",
        "Cosmic MutStat",
        "Cosmic ID"
    );

    push( @out_header_keys, qw{
            ensSIFT ensPP2hvar ensPP2hdiv ensCondel ensCondelPred
            cosmicName cosmicSite cosmicHis cosmicStat cosmicID
        } );

    # HGMD prediction
    my @hgmd_header = ( "HGMD ID" );
    push (@out_header_keys, "hgmdID");
    push (@hgmd_header, 
        "HGMD MutName",
        "HGMD Disease",
        "HGMD pmID",
        "HGMD Pred"
    );

    push( @out_header_keys, qw{
            hgmdName hgmdDis hgmdPMID hgmdPred
        } );

    # ClinVar prediction
    my @clinVar_header = ("ClinVar SuspectReason", "ClinVar Accession", "ClinVar RevStat", "ClinVar Significant");
    push (@out_header_keys, "ClinSSR", "ClinACC", "ClinRevStat", "ClinSignificant");

    # Auto Interpretation
    my @interpret_header = ("AutoInterpStatus");
    push (@out_header_keys, "autoInterp");

    my @all_headers = (
        @chrpos_header, $inExcel,           @sample_header,
        @geno_header,   @gene_header,       @func_header,
        @hgvs_header,   @freq_header,       
        @panel_header,
        @phyloP_header, @condel_header,     
        @cosmic_header, @hgmd_header,       @clinVar_header,
        @interpret_header
    );

    print $OutFp join( "\n",
        "## NCanno Version    : $VERSION",
        "## BedAnno Version   : $BedAnno::VERSION",
        "## Configure File    : " . File::Spec->rel2abs($CONFIG_FILE),
        "## Input File Format : "
          . ( ( $type ne 'tsv' ) ? "vcf" : "tsv" ),
        "## Format Options    : "
          . ( ( defined $offline )     ? 'offline, ' : 'online, ' )
          . ( ( defined $sample_list ) ? 'famlily'   : 'single' ) . " MODE",
        "## DB Version File   : " . File::Spec->rel2abs($verList),
        "##\n" );

    print $OutFp '#'.join("\t", @all_headers)."\n";

    if ($headrule) {
        print HDR '#'.join("\t", qw{title width hidden level collapsed})."\n";
        foreach my $title (@all_headers) {
            next if ($title eq "InExcel" or $title eq "FamID");
            next if (!defined $sample_list and $title eq "SampleID");
            
            # assign default stat
            my ($width, $hidden, $level, $collapsed) = (10, 1, 1, 0);

            if (   $title eq "#Chr"
                or $title eq "Ref"
                or $title eq "Call"
                or $title eq "VarType"
                or $title eq "A.Depth"
                or $title eq "A.Ratio"
                or $title eq "A.Index"
                or $title eq "Strand" )
            {
                $width = 6;
            }
            elsif ($title eq "Sex"
                or $title eq "Primary" )
            {
                $width = 6;
                ( $hidden, $level, $collapsed ) = ( 0, 0, 1 );
            }
            elsif ($title eq "Flank"
                or $title eq "RepeatTag"
                or $title eq "TestCode"
                or $title eq "FuncRegion"
                or $title eq "ExIn_ID"
                or $title eq "pfamId"
                or $title eq "EntrezGeneID" )
            {
                $width = 8;
            }
            elsif ($title eq "Filter"
                or $title eq "Zygosity"
                or $title eq "MapLoc" )
            {
                $width = 8;
                ( $hidden, $level, $collapsed ) = ( 0, 0, 1 );
            }
            elsif ($title eq "MIM Inheritance"
                or $title eq "Function"
                or $title eq "rsID"
                or $title eq "PVFD AF"
                or $title eq "Panel AlleleFreq"
                or $title eq "Panel Pred"
                or $title eq "AutoInterpStatus" )
            {
                ( $hidden, $level, $collapsed ) = ( 0, 0, 1 );
            }
            elsif ($title eq "Transcript"
                or $title eq "Protein"
                or $title eq "pHGVS" )
            {
                $width = 12;
            }
            elsif ( $title eq "Gene Symbol" ) {
                $width = 12;
                ( $hidden, $level, $collapsed ) = ( 0, 0, 1 );
            }
            elsif ( $title =~ /^HGMD |^ClinVar / ) {
                $width = 12;
                if (   $title eq "HGMD Pred"
                    or $title eq "ClinVar Significant" )
                {
                    ( $hidden, $level, $collapsed ) = ( 0, 0, 1 );
                }
            }
            elsif ($title =~ /^dbSNP |^1000G |^ESP6500 |^ExAC |^GAD /
                or $title =~ /^PhyloP |^Ens |^Cosmic / )
            {
                $width = 15;
                if (   $title eq "dbSNP Mul Tag"
                    or $title eq "1000G AF"
                    or $title eq "ESP6500 AF"
                    or $title eq "ExAC AF"
                    or $title eq "GAD AF"
                    or $title eq "PhyloP Placental Mammals"
                    or $title eq "Ens Condel Pred"
                    or $title eq "Cosmic ID" )
                {
                    $collapsed = 1;
                }
                else {
                    $level = 2;
                }
            }
            elsif ( $title eq "MutationName" ) {
                $width = 30;
                ( $hidden, $level, $collapsed ) = ( 0, 0, 1 );
            }

            print HDR join("\t", '"'.$title.'"', $width, $hidden, $level, $collapsed)."\n";
        }

        print HDR join("\t", '"Interpretation"', 15, 0, 0, 1 );
        close HDR;
    }
}

sub normalise_null {
    my @result = @_;
    @result = map {
        if ( !defined $_ ) { $_ = '.'; }
        s/^\\N$|^N\/A$|^\s*$|^null$/\./i;
        $_
    } @result;
    return @result;
}

# check the predict summary.
sub judge_pred {
    my $mPred = shift;
    $mPred =~ s/\.//g; # remove no prediction record
    my $total = length($mPred);
    my $dnum = ($mPred =~ tr/D//);
    if ($dnum >= ($total / 2.0)) {
        return "D";
    }
    else {
        return "P";
    }
}

# Reverse and complement
sub revcom {
    my $Seq = shift;
    $Seq = reverse($Seq);
    $Seq =~ tr/ATCGatcgRYKMSWBDHVrykmswbdhv/TAGCtagcYRMKSWVHDByrmkswvhdb/;
    return $Seq;
}

# read gene symbol's test code
sub read_TestCode {
    my $file = shift;
    my $TCD_h;
    if ($file =~ /\.gz$/) {
        $TCD_h = new IO::Uncompress::Gunzip $file, AUTOCLOSE=>1
            or die "Error: [$file] $GunzipError\n";
    }
    else {
        open($TCD_h, $file) or die "Error: [$file] $!\n";
    }
    my %geneTCD = ();
    while (<$TCD_h>) {
        next if (/^#|^\s*$/);
        s/\s+$//g;
        my @itm = split(/\s+/, $_, 2);
        next if (1 == @itm);
        $geneTCD{$itm[0]} = $itm[1];
    }
    close ($TCD_h);
    return \%geneTCD;
}

sub read_HomoRef_var {
    my $file = shift;
    my $HA_h;
    if ($file =~ /\.gz$/) {
        $HA_h = new IO::Uncompress::Gunzip $file, AUTOCLOSE=>1
            or die "Error: [$file] $GunzipError\n";
    }
    else {
        open($HA_h, $file) or die "Error: [$file] $!\n";
    }
    my %homoRef_var = ();
    while (<$HA_h>) {
        s/\s+$//g;
        next if (/^\s*#|^\s*$/);
        my ($chr, $start, $end, $ref) = (split(/\s+/, $_, 5))[0 .. 3];
        if (!defined $ref) {
            confess "Error format [$file]";
        }
        my $cur_var = {
            chr               => $chr,
            begin             => $start,
            end               => $end,
            referenceSequence => $ref,
            variantSequence   => $ref,
        };
        $homoRef_var{$chr} = [] if (!exists $homoRef_var{$chr});
        push (@{$homoRef_var{$chr}}, $cur_var);
    }
    close ($HA_h);
    foreach my $chrom (keys %homoRef_var) {
        my @sorted_vars = sort {
            $a->{begin} <=> $b->{begin}
         or $a->{end}   <=> $b->{end}
        } @{$homoRef_var{$chrom}};
        $homoRef_var{$chrom} = \@sorted_vars;
    }
    return \%homoRef_var;
}

# read sample infos.
sub read_sample_list {
    my $f = shift;
    my %sample_info = ();
    open (F, $f) or confess "Error: [$f] $!\n";
    while (<F>) {
        s/\s+$//g;
        next if (/^\s*#|^\s*$/);
        my @all_ids = split(/\s+/);
        $sample_info{$all_ids[0]}{sex} = $all_ids[1];
        $sample_info{$all_ids[0]}{fam} = $all_ids[2];
    }
    close F;
    return \%sample_info;
}

sub read_inmsqc {
    my $msqc_result = shift;
    my %msqc_input = ();
    open (MSQCRST, $msqc_result) or die "Error: [$msqc_result] $!\n";
    while (my $qcrst = <MSQCRST>) {
        next if ($qcrst =~ /^#/);
        chomp($qcrst);
        my @itms = split(/\s+/, $qcrst);
        my $sampid = shift( @itms );
        $msqc_input{$sampid} = [map {uc($_)} @itms];
    }
    close MSQCRST;
    return \%msqc_input;
}

sub read_MSQC_info {
    my $msqc_info_file = shift;
    my %msqc_sites = ();
    open (MSQCINFO, $msqc_info_file) or die "Error: [$msqc_info_file] $!\n";
    while (my $qcl = <MSQCINFO>) {
        next if ($qcl =~ /^#/);
        chomp($qcl);
        my @itms = split(/\t/, $qcl, 5);
        my $sid = $itms[0];
        $sid =~ s/^.+\.//;
        if ($sid =~ /^\d+$/) {
            $sid -= 1;
        }
        else {
            die join("\n", "Error: [$msqc_info] ID format in MSQC sites list",
                "should be like 'XXXX.YY', whose 'Y' should be all number\n");
        }
        $msqc_sites{$itms[1]}{$sid} = {
            pos => $itms[2],
            ref => uc($itms[3]),
        };
    }
    close MSQCINFO;
    return \%msqc_sites;
}

sub check_dep {
    foreach (@_) {
        confess "Error: [$_] no exists or can not be read or empty."
          if ( !-e $_ or !-r $_ or ( -f $_ and -z $_ ) );
    }
    return 1;
}

__END__

=head1 NAME

NCanno

=head1 VERSION

1.1.1

=head1 SYNOPSIS

    NCanno.pl <config> [options] <vcf[.gz]|tsv.bz2> ...

=head2 OPTIONS

    Optional:
        -t, --type          [STR]   type of input files vcf/tsv, default vcf,
                                    all input files should be in the same 
                                    type, mixed cases are not supported.
        -s, --samplelist    [FILE]  sample list used in pipeline, if sample
                                    list is specified, extra field "FAM_ID",
                                    "SAMPLE_ID", "Sex" will be added into 
                                    output, default only "SAMPLE_ID" got
                                    from vcf col header or from CG's input 
                                    file name is added.
        -o, --outfile       [FILE]  Output file, default STDOUT.
        -r, --headrule      [FILE]  Output header format file, default 
                                    no output.
        -i, --inmsqc        [FILE]  Result of MS Genotyping result for QC sites
        -c, --qcresult      [FILE]  Output MSQC result to FILE, default STDERR
        -n, --numthreads    [NUM]   Number of threads to be used for annotation
        -b, --buffersize    [NUM]   Number of record to be cached, ready for 
                                    sorting and output.
        -w, --writever              write DB version information to verList, 
                                    default no write.
        -f, --offline               Using an offline version of database 
                                    for emergency, network error or any other 
                                    special situation.
        -d, --debug                 Show debug information
        -m, --timer                 Show timer information
        -q, --quiet                 Suppress the warning information, default 
                                    not suppress.
        -h, --help                  Show full help message.

=head1 DESCRIPTION

    This script is used to annotate vcf file or cg's ASM-var/masterBeta tsv file
    by using BedAnno and Vcf module, and connect to various databases.

=head1 AUTHOR

Liu Tao, E<lt>liut@geneplus.org.cnE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2015 by GenePlus.

=cut

# vim: sw=4 ts=4 ft=perl expandtab
