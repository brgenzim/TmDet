<?php

require_once 'vendor/autoload.php';

use Unitmp\TmdetTest\Constants\FileSystem;

const CMD = "../build/src/tmdet -i %s -x o.xml 2>&1";

$lines = file('whole_para_2.log');
$lines = preg_grep('/non-standard/', $lines);
$lines = array_slice($lines, 0, 20);

$failedResidues = [];
$allowedResidues = [];

foreach ($lines as $line) {
    if (preg_match("/InvalidArgumentException: (\S+).+?'(\S+)'/", $line, $matches)) {
        $pdbCode = $matches[1];
        $residue = $matches[2];

        if (in_array($residue, $allowedResidues)) {
            continue;
        }

        $cifFile = FileSystem::PDB_ZFS_DIR . '/' . substr($pdbCode, 1, 2) . "/$pdbCode.cif.gz";
        $cmd = sprintf(CMD, $cifFile);
        exec($cmd, $output, $resultCode);
        $segfault = preg_grep('/core dumped/', $output) !== false;
        if ($segfault) {
            echo $cmd . PHP_EOL;
            echo $residue . PHP_EOL;
            var_dump($output);die('D');
        }
        if ($segfault && !array_key_exists($residue, $failedResidues)) {
            $failedResidues[$residue] = "gdb --args $cmd";
        } elseif (!$segfault && !in_array($residue, $allowedResidues)) {
            $allowedResidues[] = $residue;
        }
    }
}
sort($failedResidues);
sort($allowedResidues);
file_put_contents('failed_residues.ser', serialize($failedResidues));
file_put_contents('allowed_residues.ser', serialize($allowedResidues));
var_dump($allowedResidues);

