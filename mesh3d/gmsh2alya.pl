#!/usr/bin/perl

use feature "switch";

#import module for command line options
use Getopt::Long;
use Math::Trig;

$HelpDescription = " Wrong parameters.\nValid arguments are:\n 
gmsh_input_file [Alya_output_filename] [--detect] [--bcs=no|boundaries|nodes|all] [--bulkcodes]  \n 
... \n
    (older options, treat carefully) \n
    [--dontfliporder] [--nogeobounds] [--function code=FunctionFile.fun] [--dims=2|3]";

(@ARGV > 0) or die "$HelpDescription";

#Default value for options
$Dimensions='';
$TypeOfBoundaryCodes='no';
$InteriorCodes='';
$DoNotWriteBoundariesInGeo='';
$WriteBoundariesInFiles='';
%FunctionCodesAndFiles=();
$fileIn=$ARGV[0];
$fileOut=$ARGV[0];
my $DetectMode='';
$DoNotFlipNodeOrder='1';

GetOptions ("out=s" => \$fileOut, "dims=s" => \$Dimensions, "bcs=s" => \$TypeOfBoundaryCodes, 
	    "function=s" => \%FunctionCodesAndFiles, "detect"  => \$DetectMode, "nogeobounds" => \$DoNotWriteBoundariesInGeo, 
            "filebounds" => \$WriteBoundariesInFiles, "dontfliporder" => \$DoNotFlipNodeOrder, "bulkcodes" => \$InteriorCodes) 
    or die "$!, $HelpDescription";  

($TypeOfBoundaryCodes=='no' or $TypeOfBoundaryCodes=='nodes' or $TypeOfBoundaryCodes=='boundaries' or $TypeOfBoundaryCodes=='all') or die "Unrecognized option for bcs : $TypeOfBoundaryCodes";

my @NodePositions;
my @Elements0D;
my @Elements1D;
my @Elements2D;
my @Elements3D;
my %PhysicalNames;


#input file is typically something.msh
open(IN,"<$fileIn.msh") or die "I could not find input file: $fileIn.msh";

print "\n\nParsing $fileIn.msh...\n\n";
print "Detection Mode only (I will not write any output files)\n\n" if ($DetectMode);

 SECTIONS: while (my $SectionHeader = <IN>) {
     
     if ($SectionHeader =~ m/MeshFormat/) {
       MESHFORMAT: while (<IN>){
	   last MESHFORMAT if m/EndMeshFormat/;
	   chomp;
	   my @version = split;
	   ($version[0] == "2.2") or
	       die "This parser can only understand version 2.2 of the Gmsh file format";
       }
     } elsif ($SectionHeader =~ m/Nodes/) {
	 #First line should be number of nodes
	 my $NumberOfNodes = <IN>;
       NODES: while (<IN>){
	   last NODES if m/EndNodes/;
	   chomp;
	   push @NodePositions, [ split ];
       }
	 ($#NodePositions == $NumberOfNodes-1) or
	     die "Acquired a number of nodes $#NodePositions different from the one stated in the input file,  $NumberOfNodes.";
     } elsif ($SectionHeader =~ m/Elements/) {
	 #First line should be number of elements
	 if ( !defined($CombinedNumberOfElements = <IN>)) {
	     die "Found end of file when searching for elements";
	 }
       ELEMENTS: while (<IN>){
	   last ELEMENTS if /EndElements/;
	   chomp;
	   #format is
	   # elm_number  elm_type  Ntags  tag1 ... tagN node1 .. nodeL
	   my @values = split ;
	   $numberOfTags = $values[2];
	   given ($values[1]) {
	       when (15){		#Zero dimensional element:
		   #15: 1-node point.
		   push @Elements0D, [ @values[2..$#values] ];
	       }
	       when ([1, 8, 26, 27, 28] ){ #One dimensional elements:
		   #1: 2-node line.
		   #8: 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
		   #26: 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
		   #27: 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
		   #28: 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
		   push @Elements1D, [ @values[2..$#values] ];
	       }
	       when ([2,3, 9,10,16,20,21,22,23,24,25] ){ #Two dimensional elements:
		   #2: 3-node triangle.
		   #3: 4-node quadrangle.
		   #9: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
		   #10: 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
		   #16: 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
		   #20: 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
		   #21: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
		   #22: 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
		   #23: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
		   #24: 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
		   #25: 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
		   push @Elements2D, [ @values[2..$#values] ];
	       }
	       when ([4,5,6, 7, 11,12,13, 14,17,18,19,29,30,31 ]){ #Three dimensional elements:
		   #4: 4-node tetrahedron.
		   #5: 8-node hexahedron.
		   #6: 6-node prism.
		   #7: 5-node pyramid.
		   #11: 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
		   #12: 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
		   #13: 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
		   #14: 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
		   #17: 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
		   #18: 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
		   #19: 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
		   #29: 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
		   #30: 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
		   #31: 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
		   push @Elements3D, [ @values[2..$#values] ];
	       }
	       default{ die "I could not understand the element type: $values[1]";}
	   }
       }
	 ($CombinedNumberOfElements==(@Elements0D+@Elements1D+@Elements2D+@Elements3D))  or
	     die "Number of processed elements does not correspond to expected number, aborting.";
     }
     #Physical names are used for the boundaries 
     #-- CAREFUL, gmsh will only write in the mesh the nodes that belong to a physical entity.
     elsif ($SectionHeader =~ /PhysicalNames/) {
	 #First line should be number of physical entities
	 if ( !defined($NumberOfPhysicalNames = <IN>)) {
	     die "Found end of file when searching for number of physical names";
	 }
       PHYSICAL: while (<IN>){
	   last PHYSICAL if /EndPhysicalNames/;
	   chomp;
	   my @line = split;
	   $PhysicalNames{$line[1]}=$line[2];
       }
     }
     #------ NOTE: For now the parser ignores the following sections:
     elsif ($SectionHeader =~ /NodeData/) {
	 while (<IN>) {
	     last if /EndNodeData/;
	 }
     } elsif ($SectionHeader =~ /\$ElementData/) {
	 while (<IN>) {
	     last if /EndElementData/;
	 }
     } elsif ($SectionHeader =~ m/ElementNodeData/) {
	 while (<IN>) {
	     last if /EndElementNodeData/;
	 }
     } else {
	 die "Unsupported section name: $SectionHeader";
     }
}

close(IN);

############################################################
##### FINISHED READING
###########################################################


############################################################
##### WRITE OUTPUT
###########################################################

#Try to determine if it is a 2D or 3D problem
if (@Elements3D > 0) {
    print "Found some 3D elements, assuming a three dimensional geometry...\n";
    $NumberOfDimensions = 3;
    $Elements = \@Elements3D;
    $Boundaries = \@Elements2D;
} else {
    print "Could not find any 3D elements, assuming a two dimensional geometry...\n";
    print "Warning: I will drop the Z coordinate of all nodes\n\n";
    $NumberOfDimensions = 2;
    $Elements = \@Elements2D;
    $Boundaries = \@Elements1D;
}




unless ($TypeOfBoundaryCodes eq "no") {
    print "WARNING: Automated boundary code detection...\n Make sure you have ALL nodes (even non-boundary nodes) assigned to in physical entities in gmsh\n\n";

    if ($TypeOfBoundaryCodes eq "nodes" or $TypeOfBoundaryCodes eq "all") {%CodeOnNodes = ();}
    if ($TypeOfBoundaryCodes eq "boundaries" or $TypeOfBoundaryCodes eq "all") {%CodeOnBoundaries = ();}
    my %DetectedBoundaryCodes=();
    my $bound=0;
    for $row (@$Boundaries) {
	$bound++;
	my $offset = @$row[0];
	$Code = @$row[1];
	$DetectedBoundaryCodes{ $Code }=1 unless exists  $DetectedBoundaryCodes{ $Code };
	if ($TypeOfBoundaryCodes eq "nodes" or $TypeOfBoundaryCodes eq "all")
	{
	    @nodes = @$row[$offset+1..$#{$row}];
	    for my $node (@nodes) {
		$CodeOnNodes{$node}=$Code;
	    }
	}
	elsif ($TypeOfBoundaryCodes eq "boundaries" or $TypeOfBoundaryCodes eq "all")
	{
	    $CodeOnBoundaries{$bound}=$Code;
	}
    }
    $DetectedCodes=keys(%DetectedBoundaryCodes );
} else 
{
    print "Boundary codes are selected by regions, \nmake sure you have that in your Alya files...\n";
} 

print "Detecting codes in bulk elements\n";

my %DetectedBulkCodes=();

for $row (@$Elements) {
    $Code = @$row[1];
    $DetectedBulkCodes{ $Code }=1 unless exists  $DetectedBulkCodes{ $Code };
}


unless ($DetectMode) {

    open(GEO,">$fileOut.geo.dat");

    print GEO "NODES_PER_ELEMENT\n";
    $i=0;
    for $row (@$Elements) {
	$i++;
	print GEO $i," ",$#{$row}-@$row[0],"\n"; #row[0] contains the number of tags
	#   if ($#{$row}-@$row[0] != $NumberOfDimensions+1) {print "\n All hell breaks loose at $i with $#{$row}-@$row[0] nodes";}
    }
    print GEO "END_NODES_PER_ELEMENT\n";

    print GEO "ELEMENTS\n";
    $i=0;
    for $row (@$Elements) {
	$i++;
	if ($DoNotFlipNodeOrder) {
	    @nodes = @$row[@$row[0]+1..$#{$row}];
	} else {
	    @nodes = reverse(@$row[@$row[0]+1..$#{$row}]);
	}
	print GEO "$i @nodes\n"; 
    }
    print GEO "END_ELEMENTS\n";
    $DetectedElements=$i;

    print GEO "COORDINATES\n";
    $i=0;
    for $node (@NodePositions) {
	$i++;
	print GEO @$node[0]," @$node[1..$NumberOfDimensions]\n";
    }
    print GEO "END_COORDINATES\n";
    $DetectedNodes=$i;

    print GEO "BOUNDARIES\n";
    unless ($DoNotWriteBoundariesInGeo) {
	$i=0;
	for $row (@$Boundaries) {
	    $i++;
	    if ($DoNotFlipNodeOrder) {
		@nodes = @$row[@$row[0]+1..$#{$row}];
	    } else {
		@nodes = reverse(@$row[@$row[0]+1..$#{$row}]);
	    }
	    print GEO $i," @nodes\n" unless $TypeOfBoundaryCodes eq "no"; 
	}
	$DetectedBoundaries=$i;
    }
    else {
	print "\nWARNING --No boundaries will be included in the geo.dat file\n";
    }
    print GEO "END_BOUNDARIES\n";
    
    print GEO "SKEW_SYSTEMS\n";
    print GEO "END_SKEW_SYSTEMS\n";

    close(GEO);

    if ($WriteBoundariesInFiles) {
	$AccumulatedBoundaries=0;
	for my $code (sort { $a <=> $b } keys %PhysicalNames) {
	    open(BND,">$fileOut.bnd.$code.dat");
	    $i=0;
	    for $row (@$Boundaries) {
		$i++;
		if ($DoNotFlipNodeOrder) {
		    @nodes = @$row[@$row[0]+1..$#{$row}];
		} else {
		    @nodes = reverse(@$row[@$row[0]+1..$#{$row}]);
		}
		if ($CodeOnNodes{$i} eq $code)   #CAREFUL THIS IS NOT CONSISTENT WITH CODES ON BOUNDARIES
		{ 
		    print BND $i," @nodes\n";
		    $AccumulatedBoundaries++;
		}
	    }
	}
	if ($AccumulatedBoundaries ne $DetectedBoundaries) {print "WARNING: Boundaries in files do not add up to total number of boundaries\n";}
    }


    if ($InteriorCodes) {
	open(MAT,">$fileOut.mat.dat");
	$NumberOfBulkCodes = scalar( keys %DetectedBulkCodes);
	print MAT "MATERIALS, DEFAULT=1, NUMBER=$NumberOfBulkCodes\n";
	$i=0;
	for $row (@$Elements) {
	    $i++;
	    $Code =  @$row[1];
	    print MAT "$i $Code\n"; 
	}
	print MAT "END_MATERIALS\n";
    }   


    if ($TypeOfBoundaryCodes eq "nodes" or $TypeOfBoundaryCodes eq "all") {
	open(FIX,">$fileOut.fix.nod");
	print FIX "ON_NODES\n";
	for my $node (sort { $a <=> $b } keys %CodeOnNodes) {
	    print FIX "$node $CodeOnNodes{$node}\n" unless (exists $FunctionCodesAndFiles{ $CodeOnNodes{$node} } or (not exists $PhysicalNames{$CodeOnNodes{$node}}) ); 
	}
	print FIX "END_ON_NODES\n";
	close FIX;
	
	#Check for time dependent boundary values
	unless (%FunctionCodesAndFiles==()) {
	    open(TIM,">$fileOut.val.dat");
	    
	    for my $functionCode (sort { $a <=> $b } keys %FunctionCodesAndFiles) {

		#Parse the function file and store functions for each variable
		my @BoundaryFunctions=();
		open(TIMEIN,"<$FunctionCodesAndFiles{$functionCode}") or die "I could not find time function file: $FunctionCodesAndFiles{$functionCode}";
		while (<TIMEIN>) {
		    push (@BoundaryFunctions,$_);
		}
		
		my $numFunctions=@BoundaryFunctions+1;
		close(TIMEIN);
		if ( $numFunctions=!$NumberOfDimensions) {
		    die "Number of functions $numFunctions does not match number of dimensions $NumberOfDimensions.";
		}

		#Reset the hash of values in nodes 
		my %ValuesOnNodes = ();
		#Compute and store values of each boundary node that corresponds to this code
		for $row (@$Boundaries) {
		    if ($DoNotFlipNodeOrder) {
			@nodes = @$row[@$row[0]+1..$#{$row}];
		    } else {
			@nodes = reverse(@$row[@$row[0]+1..$#{$row}]);
		    }
		    for my $node (@nodes) {
			# Check if node is on the Physical domain assigned to the function
			if ($CodeOnNodes{$NodePositions[$node][0]} == $functionCode ) {
			    # Check if node has already been processed
			    unless (exists  $ValuesOnNodes{$NodePositions[$node][0] }) { 
				$X=$NodePositions[$node][1];
				$Y=$NodePositions[$node][2];
				$Z=$NodePositions[$node][3];
				my $Xval= eval $BoundaryFunctions[0];
				my $Yval= eval $BoundaryFunctions[1];
				my $Zval='';
				if ($NumberOfDimensions==3) {
				    $Zval= eval $BoundaryFunctions[2];
				    $nextIndex=3;
				} else {	  
				    $nextIndex=2;
				}
				$rhoVal= eval $BoundaryFunctions[$nextIndex];
				$TVal= eval $BoundaryFunctions[$nextIndex+1];
				$ValuesOnNodes{$NodePositions[$node][0]}="$Xval  $Yval  $Zval  $rhoVal  $TVal";
			    }
			}
		    }
		}
		# Now write the values for this function
		print TIM "VALUES, NODES, FUNCTION=$functionCode, CODE=1\n";
		for my $node (sort { $a <=> $b } keys %ValuesOnNodes) {
		    print TIM "$node $ValuesOnNodes{$node}\n"; 
		}
		print TIM "END_VALUES\n";
		
	    }
	    close(TIM);
	    
	}				#Close time dependent boundary values

    }
    if ($TypeOfBoundaryCodes eq "boundaries" or $TypeOfBoundaryCodes eq "all") {
	open(FIX,">$fileOut.fix.bou");
	print FIX "ON_BOUNDARIES\n";
	for my $bound (sort { $a <=> $b } keys %CodeOnBoundaries) {
	    print FIX "$bound $CodeOnBoundaries{$bound}\n"; 
	}
	print FIX "END_ON_BOUNDARIES\n";
	close FIX;
    }
    
    #Write velocity file:
    open(VEL,">velo.dat");
    print, "QUI>AAA";
    $i=0;
    for $node (@NodePositions) {
	$i++;
	print VEL @$node[0]," 0.0 0.0\n";
    }
    close(VEL);

    # Write dimensions file
    open (DIMS,">$fileOut.dims.dat");
    print DIMS "NODAL_POINTS    $DetectedNodes \n";
    print DIMS "ELEMENTS        $DetectedElements \n";
    print DIMS "BOUNDARIES      $DetectedBoundaries \n";
    close DIMS; 
}

print "  Statistics\n";
print "  ----------\n";

print "      Total number of nodes: $DetectedNodes\n";
print "      Total number of elements: $DetectedElements\n";
print "      Total number of boundaries: $DetectedBoundaries\n";
print "      Total number of codes for boundaries: $DetectedCodes\n";
print "         Detected codes  /   Physical Entity:\n";
print "         --------------\n";
for my $code (sort { $a <=> $b } keys %PhysicalNames) {
    print "            $code   $PhysicalNames{$code}\n";
}
unless (%FunctionCodesAndFiles==()) {
    print "\n";
    print "      The following code boundaries were assigned a value by function:\n";
    for my $functionCode (sort { $a <=> $b } keys %FunctionCodesAndFiles) {
	print "            $functionCode    $FunctionCodesAndFiles{$functionCode}\n";
    }
}


