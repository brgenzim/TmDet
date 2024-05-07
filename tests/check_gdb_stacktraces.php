<?php

$gdbCommands = unserialize(file_get_contents('failed_residues.ser'));

$checkedResidues = [];

foreach ($gdbCommands as $command) {
    $command = str_replace('gdb', 'gdb --batch -ex r -ex bt -ex exit', $command);
    exec($command, $output, $resultCode);
    $output = preg_grep('/std::unordered_map.+__k="/', $output);
    if (!$output) {
        fprintf(STDOUT, "Pattern not found in output:\n\t%s\n", $command);
    } elseif ($output && preg_match('/__k="(\S+?)"/', array_shift($output), $matches)) {
        $residue = $matches[1];
        if (!in_array($residue, $checkedResidues)) {
            $checkedResidues[$residue] = $command;
        }
    }
}


file_put_contents('gdb_checked_residues.ser', serialize($checkedResidues));
