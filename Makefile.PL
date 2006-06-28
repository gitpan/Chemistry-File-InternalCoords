use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'Chemistry::File::InternalCoords',
    AUTHOR              => 'David Westbrook <dwestbrook@gmail.com>',
    VERSION_FROM        => 'lib/Chemistry/File/InternalCoords.pm',
    ABSTRACT_FROM       => 'lib/Chemistry/File/InternalCoords.pm',
    PL_FILES            => {},
    PREREQ_PM => {
        'Test::More' => 0,
	'Chemistry::Mol' => 0,
	'Chemistry::File' => 0,
	'Chemistry::InternalCoords' => 0,
	'Chemistry::InternalCoords::Builder' => 0,
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Chemistry-File-InternalCoords-*' },
);
