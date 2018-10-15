

Using capC-MAP
==============

The capC-MAP software is essentially a collection of separate tools
along with a "wrapper script", and can be used in two ways. Either, the
full Capture-C analysis can be performed by using the single command
``capC-MAP run``, or for advanced use individual steps of the pipeline can
be run separately.

The ``capC-MAP`` command can be run in one of five mode:

   * ``genomedigest`` is used to generate a list of restriction enzyme fragments covering the reference genome; 
   * ``run`` is used to run the full analysis pipeline;
   * ``postprocess`` is used to generate binned contact profiles;
   * ``combinereps`` is used to combine processed data sets from replicates; and
   * ``getchromsizes`` is used to generate a list of chromosome sizes from a list or restriction enzyme fragments.

The ``run``, ``postprocess`` and ``combinereps`` modes read from a "configuration file"; the other modes take options on the command line.  

.. _secgendig:

capC-MAP ``genomedigest``
*************************

In order to identify interactions between restriction enzyme fragments,
capC-MAP needs a genome-wide map of where the restriction sites are. Using
the ``genomedigest`` mode capC-MAP can generate
this list from a fasta file for the reference genome. Once generated
this can be re-used for any Capture-C experiments using that reference genome.

capC-MAP stores the list of fragments as a standard bed file; this can
be generated with the following command:

::

   capC-MAP genomedigest -f <reference genome fasta file> \
                         -r <name of restriction enzyme|cut sequence> \
                         -o <output bed file>

where all three options are required. The ``-r`` option specifies either the
name of a restriction enzyme which is known to capC-MAP (i.e. DpnII), or the sequence recognised by the enzyme. Typically DpnII is used in Capture-C experiments, but the ``-r`` option allows other cut sequences to be used (see section :ref:`sectionenzyme` below). 

.. _sectionrun:

capC-MAP ``run``
****************

The capC-MAP software performs each of the following processing steps
using a single command line:

* trim adapters from paired-end read fastq files using the **cutadapt** software [Martin2011]_
* perform *in silico* restriction enzyme digestion,
* align restriction enzyme fragments to reference genome using **bowtie** [Langmead2009]_ in single-end read mode,
* sort the resulting SAM file by read name using **samtools** [Li2009]_,
* identify “read groups” of mapped restriction enzyme fragments and remove duplicates,
* identify “target” and “reporter” fragments in each group,
* remove interactions between targets, interactions within “exclusion zones”, and groups with multiple non-adjacent reporters,
* remove and count interchromosomal interactions,
* generate pile-ups of interactions for each target,
* generate normalized, binned and smoothed interaction profiles for each target.

Full details of each of these steps are given in the associated
paper [capC-MAP]_. This pipeline is based on that
detailed in the original Capture-C publications from the Hughes
Lab [Hughes2014]_ (and see also [Davies2016]_).

The command line takes the form:

::

  capC-MAP run -c <configuration file> \
               -o <output directory>   

where all options and input files are specified in the configuration file, and capC-MAP will create a new directory at the
specified location for all output files. capC-MAP will not overwrite existing files, and will instead fail with an error.
	       
The configuration file
----------------------

The configuration file is a text file which sets all of the
options required for the capC-MAP pipeline. We recommend taking the
example configuration file provided with the capC-MAP installation as a
template. The file is structured such that each line represents an
option specified by the first word on that line. Any line beginning with
a hash ‘#’ character is treated as a comment and ignored by capC-MAP.
Some options are required to be in the configuration file, whereas others
will take a default value if not specified.

The available options are as follows

``FASTQ1 <file>``
  *Required*. Specifies the relative path to the first of the pair of fastq files. Takes exactly one argument; subsequent arguments are ignored.

``FASTQ2 <file>``
  *Required*. Specifies the relative path to the second of the pair of fastq files. Takes exactly one argument; subsequent arguments are ignored.

``TARGETS <file>``
  *Required*. Specifies a bed file containing a list of “target”
  restriction enzyme fragments. Each target fragment must match as
  fragment given in the ``RESTFRAGS`` file. Details are given in
  section :ref:'sectiontargets' below. Takes exactly one argument;
  subsequent arguments are ignored.

``INDEX <index>``
  Required. Specifies the bowtie index for the reference genome;
  this will be passed to bowtie for the alignment step of the
  pipeline. The relative path must be given, with the name as
  specified when the reference was built using ``bowtie-build``.
  Takes exactly one argument; subsequent arguments are ignored.

``RESTFRAGS <file>``
  *Required*. Specifies a bed file for the list of restriction enzyme
  fragments covering the reference genome, as generated using
  capC-MAP with the ``genomedigest`` mode. Note that this must
  have been generated for the same restriction enzyme as is specified
  by the ``ENZYME`` option. Takes exactly one argument; subsequent
  arguments are ignored.

``ENZYME <enzyme name|cut sequence>``
  *Required*. Specifies the restriction enzyme used in the experiment, either by name or by 
  specifying the recognition sequence directly.
  As detailed in section :ref:`secgendig` above, this is typically DpnII in a Capture-C 
  experiment, but alternatively the recognition sequence can
  be specified directly. See section :ref:`sectionenzyme` below for further details.
  Takes exactly one argument, case insensitive; subsequent
  arguments are ignored.

``TRIMADAPTERS [TRUE|FALSE]``
  *Optional*. Default: TRUE. Since the Capture-C protocol recommends
  sonicating the library to give short fragments, it is expected that
  for many reads sequencing will have gone through into the adapter
  sequence. By default capC-MAP uses the cutadapt software to trim
  adapters from the input fastq files. This step can be skipped by
  setting this option to FALSE. Takes exactly one argument; subsequent
  arguments are ignored.

``PARALLEL <N>``
  *Optional*. Default: 1. To speed up processing, some step of
  the capC-MAP pipeline can be run on multiple processors. Specifically 
  sequence alignment using bowtie and sorting and file conversion using 
  samtools can be run in parallel. This option specifies the number of 
  processors, and this is passed to bowtie and samtools. Since these are the 
  slowest steps in the pipeline, no other steps are run on multiple processors. 
  Takes exactly one integer argument; subsequent arguments are ignored.

``ALIGNMODE [CONSERVATIVE|RELAXED|CUSTOM]``
  *Optional*. Default: CONSERVATIVE. Determines the alignment options
  which are passed to bowtie. See section :ref:`secalignmode`
  below for details. The two pre-set options CONSERVATIVE and RELAXED
  require no further options. If CUSTOM is specified, everything
  following it on the same line is taken to be an option for the
  aligner, and is passed verbatim to bowtie.

``EXCLUDE <N>``
  *Optional.* Default: 500. Sets the distance in base-pairs, where if a 
  reporter fragment is closer to a target than this it is discarded.

``INTERCHROM``
  *Optional.* Default: FALSE. Sets whether capC-MAP generates pile-ups for 
  interchromosomal interactions. Note that valid interchromosomal interaction 
  pairs are always saved in 'validinterchrom' output files.

``BIN <S> <W>``
  *Optional.* Tells capC-MAP to generate binned interaction profiles as well 
  as restriction enzyme fragment level pile-ups. Since restriction enzyme 
  fragments have an irregular size, some binning is recommended. The step size 
  for bins in base-pairs is set by the integer ``<S>``. Smoothing can also be 
  applied via a sliding window of width ``<W>`` base-pairs, i.e. each bin 
  gives the number of interactions from within a window of that width. To bin 
  without smoothing set W=S. If the ``NORMALIZE`` option is also set TRUE, the 
  binned profiles will be normalized.

``NORMALIZE``
  *Optional.* Default: FALSE. Sets whether binned interaction profiles are 
  also normalized to be in units of "reads-per-million".

``COMBINEMODE [TRUE|FALSE]``
  *Optional*. Default: FALSE. Specifies whether the analysis will be
  run with ‘combine mode’ activated. This will combine interactions
  from selected targets into a single interaction profile. See
  section :ref:`seccombine` below for details. Takes exactly
  one argument; subsequent arguments are ignored.

``COMBINECOUNT <N>``
  *Optional*. Default: 2. Only relevant when combine mode is active.
  Specifies how many targets are to be combined. See
  section :ref:`seccombine` below for details. Takes exactly
  one integer argument; subsequent arguments are ignored.

``DRYRUN [TRUE|FALSE]``
  *Optional*. Default: FALSE. If set TRUE capC-MAP will be run in "dry run"
  mode, which steps through each stage of the pipe-line without actually
  running it. This is useful for testing all required files etc. are present,
  and generating the 'capC-MAP.commands.log' file, which lists all pipe-line
  steps as bash command lines (see section :ref:`secoutputs` below).

.. _secalignmode:

Alignment mode
--------------

Alignment of fragments to the reference genome is done using the bowtie
software. By default capC-MAP uses a rather conservative set of
parameters for the alignment, as recommended in the original Capture-C
protocol [Hughes2014]_: only fragments which map to a
single genomic location are reported. Another, more relaxed, pre-set for
the alignment parameters where the best alignment for multi-mapping reads
is reported is also available. Alternatively users can specify their own
custom set of parameters to be passed to bowtie. Note that since the
number of processors which bowtie uses is specified separately, this
should not be included in the custom alignment mode line of the
configuration file.

.. _sectiontargets:

The targets file
----------------

A bed file containing a list of all targeted restriction enzyme fragments is a required input, and is specified with the ``TARGETS`` option in the configuration file. Each line must contain four fields separated by tabs: chromosome, start, end, and target name; for example

::
   
   chr2     12345     67890     firsttargetname
   chr4     23456     78901     secondtargetname

Target names must be unique, and the same fragment cannot appear more than once. Each target fragment must be also present in the restriction enzymes file specified with the ``RESTFRAGS`` option in the configuration file (i.e. there must be a line with the same chrom:start:end fields). capC-MAP provides a utility ``location2fragment`` which is useful for generating a valid targets file, e.g. from a bed file containing a list of oligo regions - see section :ref:`sectionloc2frag`.

.. _sectionenzyme:

Restriction Enzymes
-------------------

capC-MAP needs to know the restriction enzyme used to digest the genome. For Capture-C experiments this is typically DpnII. To use a different restriction enzyme, the DNA recognition sequence can be specified with the ``ENZYME`` option. Note that for recognition sequences with an overhang, such as HindIII, there is a nucleotide fill-in step during ligation - this means that the recognition sequence used to generate the restriction enzyme map with the capC-MAP ``genomedigest`` command may be different to the one used for the rest of the analysis.

.. _seccombine:

Combine mode
------------

Sometimes a genome feature of interest might appear at multiple locations
in the genome. For example, in the paper in which the Capture-C method
was originally described, the authors studied interactions with
the promoter of the mouse :math:`\alpha`-globin gene. There are two
copies of :math:`\alpha`-globin the mouse genome, with largely the same
sequence. While oligos designed to target those promoter will lead to
enrichment of fragments containing either copy, these will be associated
to only a single genomic location when aligned to the reference genome.
Thus interactions for the two targets should be combined into a single
interaction profile. This is handled automatically when capC-MAP is run
in “combine mode”, provided that the targets to be combined are named in
a specific way. Names of targets which start with the same string, and
end with “``_C1``”, “``_C2``”, “``_C3``” :math:`\ldots` etc. For example,
an experiment targeting the two copies of the mouse :math:`\alpha`-globin gene (mm9) might use the following targets:

::

   chr11    32182970     32183819     AGLOB_C1
   chr11    32195805     32196636     AGLOB_C2

When run with combine mode set TRUE, capC-MAP will generate a set of output 
files with "AGLOB_combined" as the target name, as well as output for 
"AGLOB_C1" and "AGLOB_C2" individually.
When capC-MAP is run with “conservative” alignment mode (recommended),
options are passed to bowtie which specify that only reads which map
uniquely to a single location are reported. When combine mode is used,
target fragments are likely to map to multiple locations, so we must
relax this restriction. This is done with the ``COMBINECOUNT`` option:
if two targets are to be combined, this should be set to 2; if three
targets are to be combined, this should be set to 3, etc. By default
``COMBINECOUNT`` will be set to 2 when combine mode is active.

.. _secoutputs:

Outputs
-------

capC-MAP generates the following files in the output directory:

capC-MAP.commands.log
  A log file showing a list of command lines for each step of the analysis.
  This is also generated in DRYRUN mode.

captured_report.dat
  A report file from the main processing stage of the pipe-line.
  Shows counts of various points where reads were discarded,
  useful for evaluating the quality of the data.
  
captured_interactioncounts.dat
  Contains counts for each target of the number of valid interactions,
  and how many were intra/inter chromosomal.

srt_aligned.bam
  BAM file for the aligned read fragments sorted by name

captured_validpairs\_\ *targetname*.pairs
  A set of files containing a list of all valid intrachromosomal interactions,
  one file for each target. Restriction enzyme fragment coordinates are given
  in bed file format.
  
captured_validinterchom\_\ *targetname*.pairs
  Similar files showing interchromosomal interactions. 
  
captured_rawpileup\_\ *targetname*.bdg
  Set of bedGraph files, one for each target, giving the "piled-up" intrachromosomal
  interactions. Each entry refers to a single restriction enzyme fragment, so these
  have irregular widths. Units are numbers of reads.
  
captured_normalizedpileup\_\ *targetname*.bdg
  When the ``NORMALIZE`` parameter is set TRUE, capC-MAP also generates a set 
  of bedGraph files where the piled-up intrachromosomal interaction counts 
  have been normalized to reads-per-million, i.e. the number of reads for 
  each target genome wide will sum to one million.

captured_rawpileup_interchom\_\ *targetname*.bdg
  Set of bedGraph files, giving the "piled-up" interchromosomal
  interactions. Only present if option
  ``INTERCHROM`` was set TRUE in the configuration file.

captured_normalizedpileup_interchom\_\ *targetname*.bdg
  When the ``NORMALIZE`` and ``INTERCHROM`` options are both set TRUE, capC-MAP also generates a set 
  of bedGraph files where the piled-up interchromosomal interaction counts 
  have been normalized to reads-per-million.

captured_bin\_\ *S*\ \_\ *W*\ \_\ *targetname*.bdg
  Here *S* ans *W* are integers. Set of bedGraph files containing the intrachromosomal
  interaction profile which has been binned using a step size *S* and a window size
  *W*, one file for each target. Units are numbers of reads.

captured_bin\_\ *S*\ \_\ *W*\ _RPM\_\ *targetname*.bdg
  As above, but units are in reads-per-million (RPM). These are generated
  instead of the above if option ``NORMALIZE`` was set TRUE in the configuration file.

Additionally log files and error files are generated from each step of the pipe-line, and these contain any output from the programs used in each step - this is useful for troubleshooting if capC-MAP fails with an error.

.. _sectionpostprocess:

capC-MAP ``postprocess``
************************

The postprocess mode is used to generate additional interaction profiles from a data set which has
already been analysed using the ``capC-MAP run`` command. By adding new ``BIN`` or ``NORMALIZE`` lines
to the configuration file, this mode can be used to generate new 
captured_bin\_\ *S*\ \_\ *W*\ \_\ *targetname*.bdg or
captured_bin\_\ *S*\ \_\ *W*\ _RPM\_\ *targetname*.bdg files from the 
captured_rawpileup\_\ *targetname*.bdg or 
captured_normalizedpileup\_\ *targetname*.bdg files.

The command line takes the form

::

  capC-MAP postprocess -c <configuration file> \
                       -o <output directory>

where the configuration file and output directory are the same ones used in the original ``capC-MAP run`` command.

If the pile-up files are not present in the output directory, capC-MAP will try to generate them from the 
captured_validpairs\_\ *targetname*.pairs files. By default capC-MAP will not overwrite any existing files.
		       


capC-MAP ``combinereps``
**************************

capC-MAP provides a facility for combining replicate data sets. Each set must first be analysed independently using the ``capC-MAP run`` command; then the ``capC-MAP combinereps`` command is used to combine the data into a single set of files for each target, and generate binned, smoothed and normalized interaction profiles. A typical set of commands might take the form

::

   capC-MAP run -c config_rep1.txt -o output_rep1
   capC-MAP run -c config_rep2.txt -o output_rep2
   capC-MAP combinereps -c config_rep1.txt  \
                        -i output_rep1 \
			-i output_rep2 \
			-o output_combinedreps_1_2

where the first two commands run the analysis on each of the replicates, and the third combines that data. The two replicate data sets must be generated using the same targets file, genome index and restriction enzyme, as specified in the configuration file. Multiple instances of the ``-i`` option are used to select the directories containing the capC-MAP output for each individual replicate. Once the combined results directory has been generated, new Capture-C profile files with different binning smoothing and normalization options can be generated using the ``capC-MAP postprocess`` command as detailed above.

capC-MAP ``getchromsizes``
**************************

It is often useful to have a list of the chromosome sizes for a reference genome, and capC-MAP provides a tool to generate this from a restriction enzyme map, as generated using the ``capC-MAP genomedigest`` command. The command line takes the form

::

   capC-MAP getchromsizes -f <fragments file> \
                          [-o <output file>]

where the fragments file must be in the format generated by ``capC-MAP genomedigest``, and if the optional output file is not specified, the file name ``chrom.sizes`` will be used.
