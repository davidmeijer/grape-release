#!/usr/local/bin/perl
open( MYFILE, '<', '../monomer_types/amino_acids.txt' );
open( OUTFILE, '>', '../monomer_types/amino_acid_domains_output.txt' );
while(<MYFILE>) {
	chomp;
	if( (m/[Aa]lanine/) && (not m/[Pp]henylalanine/) && (not m/[Bb]eta-[Aa]lanine/) ) { s/$/\tAla/g; }
	if( m/[Bb]eta-[Aa]lanine/ ) { s/$/\tbAla/g; }
	if( m/[Pp]henylalanin/ && m/bOH/ ) { s/$/\tOHPhe/g; }
	elsif( m/[Pp]henylalanin/ ) { s/$/\tPhe/g; }
	if( m/hydroxyleucine/ ) { s/$/\tOHLeu/g; }
	elsif( m/(?<!([Ii]so))[Ll]eucin/ ) { s/$/\tLeu/g; }
	if( (m/Leu/) && (not m/[Ll]eucin/) ) { s/$/\tLeu/g; }
	if( m/[Ii]soleucin/ ) { s/$/\tIle/g; }
	if( m/[Dd]i[Hh]ydroxy[Pp]henyl[Gg]lycin/ ) { s/$/\tDhpg/g; }
	elsif( m/[Hh]ydroxy[Pp]henyl[Gg]lycin/ ) { s/$/\tHpg/g; }
	elsif( m/[Pp]henyl[Gg]lycin/ ) { s/$/\tHpg/g; }
	elsif( m/[Gg]lycin/ ) { s/$/\tGly/g; }
	if( m/[Hh]ydroxy[Aa]sparagine/ ) { s/$/\tOHAsn/g; }
	elsif( m/[Aa]sparagine/ ) { s/$/\tAsn/g; }
	if( m/[Hh]ydroxy[Aa]spartic [Aa]cid/ ) { s/$/\tOHAsp/g; }
	elsif( m/[Aa]spartic [Aa]cid/ && m/bMe/) { s/$/\tMeAsp/g; }
	elsif( m/[Aa]spartic [Aa]cid/ ) { s/$/\tAsp/g; }
	if( m/2,3-dihydroxybenzoic acid/ ) { s/$/\tSal/g; }
	elsif( m/[Bb]enzoic [Aa]cid/ ) { s/$/\tBz/g; }
	if( m/[Cc]apreomycidine/ ) { s/$/\tCap/g; }
	if( m/[Cc]itrulline/ ) { s/$/\tCit/g; }
	if( m/[Cc]ysteine/ ) { s/$/\tCys/g; }
	if( m/[Cc]ysteic/ ) { s/$/\tCys/g; }
	if( m/[Gg]lutamine/ ) { s/$/\tGln/g; }
	if( m/[Gg]lutamic [Aa]cid/ ) { s/$/\tGlu/g; }
	if( m/[Hh]istidine/ ) { s/$/\tHis/g; }
	if( m/[Mm]ethionine/ ) { s/$/\tMet/g; }
	if( m/[Nn]orvalin/ ) { s/$/\tNva/g; }
	elsif( m/[Ii]sovalin/ ) { s/$/\tIva/g; }
	elsif( m/[Vv]alin/ && m/bOH/) { s/$/\tOHVal/g; }
	elsif( m/[Vv]alin/ ) { s/$/\tVal/g; }
	if( m/[Oo]rnithine/ && m/[Hh]ydroxy/ ) { s/$/\tOHOrn/g; }
	elsif( m/[Oo]rnithine/ ) { s/$/\tOrn/g; }
	if( m/[Tt]ryptophan/ ) { s/$/\tTrp/g; }
	if( m/Bmt/ ) { s/$/\tBmt/g; }
	elsif( m/[Tt]hreonine/ ) { s/$/\tThr/g; }
	if( m/[Aa]rgin/ ) { s/$/\tArg/g; }
	if( m/[Mm]ethyl[Pp]roline/) { s/$/\tMePro/g; }
	elsif( m/[Hh]omo[Pp]roline/ ) { s/$/\tPip/g; }
	elsif( m/[Pp]roline/ ) { s/$/\tPro/g; }
	if( m/[Ss]erin/ ) { s/$/\tSer/g; }
	if( m/[Tt]yrosine/ && m/[Hh]ydroxy/ ) { s/$/\tOHTyr/g; }
	elsif( m/[Tt]yrosine/ ) { s/$/\tTyr/g; }
	if( m/bLys/ ) { s/$/\tbLys/g; }
	elsif( m/[Ll]ysine/ ) { s/$/\tLys/g; }
	if( m/[Aa]dipic/ ) { s/$/\tAad/g; }
	if( m/[Ll]actate/ ) { s/$/\tLac/g; }
	if( m/[Ll]actic acid/ ) { s/$/\tLac/g; }
	if( m/[Pp]ropionyl/ ) { s/$/\tHap/g; }
	if( m/hydroxyisovalerate/ ) { s/$/\tHiv/g; }
	if( m/kynurenine/ ) { s/$/\tKyn/g; }
	if( m/[Aa]bu/ && m/dehydro/ ) { s/$/\tDhab/g; }
	elsif( m/[Aa]bu/ ) { s/$/\tAbu/g; }
	if( m/[Aa]minoisobutyric/ ) { s/$/\tAib/g; }
	if( m/[Cc]oronamic/ ) { s/$/\tCma/g; }
	if( m/[Dd]olaproine/ ) { s/$/\tDap/g; }
	if( m/[Ee]nduracididine/ ) { s/$/\tEnd/g; }
	if( m/Hmp/ ) { s/$/\tHmp/g; }
	if( m/Dpr/ ) { s/$/\tDpr/g; }
	if( m/[Vv]aleric/ ) { s/$/\tVaa/g; }
	if( m/Dab/ ) { s/$/\tDab/g; }	
	if( m/Dbu/ ) { s/$/\tDab/g; }	# diaminobutyric acid
	if( m/MethylGlutamate/ ) {s/$/\tMeGlue/g}
	if( m/[Pp]ipecolic/ ) {s/$/\tPip/g}
	if( m/[Pp]yruvate/ ) {s/$/\tPyr/g}
	if( m/isocaproate/ ) {s/$/\takIL/g}
	if( m/isovalerate/ && m/keto/ ) {s/$/\takIV/g}
	if( m/oxopentoate/ ) {s/$/\tOxo/g}

	# modifications. Count the number of instances of \t - if more than one, then go over these possible additional modifications
	my $count = 0;
	for my $c (split('')) {$count++ if $c eq "\t";}

	if($count > 1) {
#		if( m/D-/ ) { s/$/\tEPIMERIZATION/g; }
		if( m/NMe/ ) { s/$/\tNMT/g; }
		if( m/OMe/ ) { s/$/\tOMT/g; }
#		if( m/NFo/ ) { s/$/\tNFORMYL/g; }
#		if( m/Ac-/ ) { s/$/\tNACETYL/g; }
#		if( m/[Mm]ethoxy/ ) { s/$/\tMETHOXY/g; }
#		if( m/[Hh]ydroxy/ ) { s/$/\tHYDROXY/g; }
#		if( m/[Pp]h-/ ) { s/$/\tPHENYL/g; }
	}
	print OUTFILE $_ . "\n";
	}

close MYFILE;
close OUTFILE;
