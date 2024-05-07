<?php

require_once 'vendor/autoload.php';

use Unitmp\TmdetTest\Constants\FileSystem;

const CMD = "../build/src/tmdet -i %s -x o.xml 2>&1";

$lines = file('residues/value_error_codes.txt');
$failedResidues = [];
foreach ($lines as $pdbCode) {

    $pdbCode = trim($pdbCode);
    $cifFile = FileSystem::PDB_ZFS_DIR . '/' . substr($pdbCode, 1, 2) . "/$pdbCode.cif.gz";
    $cmd = sprintf(CMD, $cifFile);
    exec($cmd, $output, $resultCode);
    $segfault = preg_grep('/core dumped/', $output) !== false;
    if ($segfault) {
        $failedResidues[] = "gdb --batch -ex r -ex bt -ex q --args $cmd";
        file_put_contents('failed_residues.ser', serialize($failedResidues));
    }
}
