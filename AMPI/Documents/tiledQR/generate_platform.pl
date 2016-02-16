#!/usr/bin/perl

if($#ARGV<1){
    printf("Usage: %s PROCESSOR_COUNT PROCESSOR_SPEED > platform.xml \nPROCESSOR_SPEED denotes the FLOP/s of a single processor\n",$0);
    exit;
}


my $proc_cnt = shift;
my $proc_speed = shift;

print "<?xml version='1.0'?>\n";
print "<!DOCTYPE platform SYSTEM \"simgrid.dtd\">\n";
print "<platform version=\"3\">\n";
print "  <AS id=\"AS0\" routing=\"Full\">\n";
print "    <host id=\"master\" power=\"10000000\"/>\n";

for(my $i=0;$i<$proc_cnt;$i++){
    print "    <host id=\"p$i\" power=\"$proc_speed\"/>\n";
}

print "    <link id=\"backbone\" bandwidth=\"100000000\" latency=\"0\"/>\n";
print "    <link id=\"loopback\" bandwidth=\"498000000\" latency=\"0.0000\"/>\n";

print "\n    <route src=\"master\" dst=\"master\"><link_ctn id=\"loopback\"/></route>\n";
for(my $i=0;$i<$proc_cnt;$i++){
    print "    <route src=\"master\" dst=\"p$i\"><link_ctn id=\"backbone\"/></route>\n";
    print "    <route dst=\"master\" src=\"p$i\"><link_ctn id=\"backbone\"/></route>\n";
}
for(my $i=0;$i<$proc_cnt;$i++){
    print "\n";
    print "    <route src=\"p$i\" dst=\"p$i\"><link_ctn id=\"loopback\"/></route>\n";
#    for(my $j=0;$j<$proc_cnt;$j++){
#        if($i==$j){
#            print "    <route src=\"p$i\" dst=\"p$i\"><link_ctn id=\"loopback\"/></route>\n";
#
#        }
#        else{
#            print "    <route src=\"p$i\" dst=\"p$j\"><link_ctn id=\"backbone\"/></route>\n";
#            print "    <route dst=\"p$i\" src=\"p$j\"><link_ctn id=\"backbone\"/></route>\n";
#        }
#    }
}
print "  </AS>\n";
print "</platform>\n";
