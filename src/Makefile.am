AUTOMAKE_OPTIONS = foreign

BUILD_DIR	=	build

bin_PROGRAMS	=	$(top_builddir)/${BUILD_DIR}/capCmain	\
			$(top_builddir)/${BUILD_DIR}/capCdigestfastq	\
			$(top_builddir)/${BUILD_DIR}/capCpair2bg	\
			$(top_builddir)/${BUILD_DIR}/capCpileup2binned	\
			$(top_builddir)/${BUILD_DIR}/capClocation2fragment

__top_builddir____BUILD_DIR__capCmain_SOURCES = main_process.cc	\
				bedfiles.cc	\
				genome.cc	\
				messages.cc	\
				parse_sam.cc	\
				samfragments.cc	\
				targets.cc
__top_builddir____BUILD_DIR__capCdigestfastq_SOURCES = fqdigest.cc	\
				fastq.cc		\
				messages.cc
__top_builddir____BUILD_DIR__capCpair2bg_SOURCES = pair2bg.cc	\
				bedfiles.cc
__top_builddir____BUILD_DIR__capCpileup2binned_SOURCES = pileup2binned.cc\
					binprofile.cc	\
					bedgraphfiles.cc\
					bedfiles.cc
__top_builddir____BUILD_DIR__capClocation2fragment_SOURCES = location2fragment.cc\
					bedfiles.cc	\
					genome.cc	\
					messages.cc	\
					targets.cc

AM_CPPFLAGS = -c -O3 -Wall
