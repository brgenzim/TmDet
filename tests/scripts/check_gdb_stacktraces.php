<?php

if ($argc != 2) {
    die('Usage: php ' . $argv[0] . " stacktrace_file\n");
}

$traceLines = explode("\n", file_get_contents($argv[1]));
$traceLines = preg_grep('/std::unordered_map.+__k="/', $traceLines);


$nonStdResidues = [];

foreach ($traceLines as $line) {
    if (preg_match('/__k="(\S+?)"/', $line, $matches)) {
        $residue = $matches[1];
        if (!in_array($residue, $nonStdResidues)) {
            $nonStdResidues[] = $residue;
        }
    }
}

sort($nonStdResidues);

foreach (array_chunk($nonStdResidues, 20) as $chunk) {
    foreach ($chunk as $residue) {
        printf("%5s, ", "'" . trim($residue) . "'");
    }
    printf("\n");
}
