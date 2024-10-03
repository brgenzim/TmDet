<?php

use Unitmp\TmdetTest\PdbtmDifferencesTest;
use Unitmp\TmdetTest\Services\PdbtmComparator;
use Unitmp\TmdetTest\Services\Tmdet30Runner;

require_once 'vendor/autoload.php';

const COMPARISION_RESULT_FILE = './tmdet_differences.json';

$codes = PdbtmDifferencesTest::getEntPathsOfTmps();
$compResults = [];
foreach ($codes as $file => $item) {
    echo "Processing data of $file" . PHP_EOL;
    $runner = Tmdet30Runner::createRunner($file);
    // $runner->enableOverwrite = true;
    try {
        $runner->exec();
    } catch (Exception $e) {
        fprintf(STDERR, "EXCEPTION: %s - %s\n", $runner->pdbCode, $e->getMessage());
        continue;
    }
    echo "meld {$runner->oldTmdetFile} {$runner->newTmdetFile}" . PHP_EOL;
    $differences = PdbtmComparator::compareTmdetData($runner);
    if (empty($differences['messages'])) {
        continue;
    }
    $compResults[] = [
        'code' => $runner->pdbCode,
        'details' => $differences
    ];
}

echo 'Writing result into ' . COMPARISION_RESULT_FILE . ' ... ';
fflush(STDOUT);
file_put_contents(COMPARISION_RESULT_FILE,
    json_encode($compResults, JSON_PRETTY_PRINT|JSON_UNESCAPED_SLASHES));
echo 'DONE.' . PHP_EOL;
