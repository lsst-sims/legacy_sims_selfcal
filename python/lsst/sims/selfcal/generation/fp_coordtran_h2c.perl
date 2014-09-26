#!/usr/bin/perl -w
# script to generate a C code file given the #defines written in the coordtran
# .h files
# expect input on stdin

$input=$ARGV[0] if defined($ARGV[0]);;
while (<>) {
    if (/FP model/) {
	($modelstr)=($_ =~ /FP model (\C+)/);
	chomp $modelstr;
    }
    next if (!/define/);
    @_=split;
    push(@DEFS,$_[1]);
}

printf STDERR "number of defines: %d\n",$#DEFS;

$list=join(",\n",@DEFS);
$nlist=$#DEFS-$[+1;

# write out a .c file intended to be #included by another .c file.
# structure definition will be included in a .h file and included in the
# .c file that initializes the structure array parsing the FP model properly.

open(HDRF,">ampcoordtran.h") || die;

printf HDRF 
    "typedef struct amp_coordtran {\n".
    "int   raft_ix[2];\n".
    "int   sens_ix[2];\n".
    "int   device_ix[2];\n".
    "float xo[2];\n".
    "float xc[2][2];\n".
    "float po[2];\n".
    "float pc[2][2];\n".
    "int   addr_lim[2][2];\n".
    "float addr_elim[2][2];\n".
    "} amp_coordtran;\n";
printf HDRF "int init_amp_coordtran(amp_coordtran *acs_array,int n_acs);\n";
printf HDRF "#define FINGERPRINT \"%s\"\n",$modelstr;
printf HDRF "#define N_AMP_COORDTRAN %d\n",$nlist;
close(HDRF);

open(SRCF,">init_ampcoordtran.c") || die;

if (defined($input)) {
    printf SRCF "#include \"$input\"\n" ;
} else {
    printf SRCF "// uncomment and edit the following line appropriately:\n// #include \"coordtran_geometryspec.h\"\n";
}

printf SRCF "
#include <stdio.h>
#include \"ampcoordtran.h\"
";
printf SRCF "char *amp_coordtran_str[]={\n%s\n};\n",$list;
printf SRCF "
int init_amp_coordtran(amp_coordtran *acs_array,int n_acs) {
  int i;
  amp_coordtran *acs;
  for (i=0;i<n_acs;i++) {
    acs=&acs_array[i];
    sscanf(amp_coordtran_str[i],
           \"R%%1d%%1dS%%1d%%1dC%%1d%%1d \"
           \"%%f %%f %%f %%f %%f %%f \"
           \"%%f %%f %%f %%f %%f %%f \"
           \"%%d %%d %%d %%d \"
           \"%%f %%f %%f %%f\",
           &acs->raft_ix[0],       &acs->raft_ix[1],
           &acs->sens_ix[0],       &acs->sens_ix[1],
           &acs->device_ix[0],     &acs->device_ix[1],
           &acs->xo[0],            &acs->xo[1],
           &acs->xc[0][0],         &acs->xc[0][1],
           &acs->xc[1][0],         &acs->xc[1][1],
           &acs->po[0],            &acs->po[1],
           &acs->pc[0][0],         &acs->pc[0][1],
           &acs->pc[1][0],         &acs->pc[1][1],
           &acs->addr_lim[0][0],   &acs->addr_lim[0][1],           
           &acs->addr_lim[1][0],   &acs->addr_lim[1][1],           
           &acs->addr_elim[0][0],  &acs->addr_elim[0][1],          
           &acs->addr_elim[1][0],  &acs->addr_elim[1][1]
           );
  }
  return(0);
}
";
close(SRCF);
