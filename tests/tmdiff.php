<?php

use Unitmp\TmdetTest\PdbtmDifferencesTest;
use Unitmp\TmdetTest\Services\Tmdet30Runner;

require_once 'vendor/autoload.php';

$codes = PdbtmDifferencesTest::getEntPathsOfTmps();
foreach ($codes as $file => $item) {
    echo $file . PHP_EOL;
    // if (substr($file, 3, 4) !== '1ldf') {
    //     printf("{$item[0]}\n");
    //     continue;
    // }
    // $runner = Tmdet30Runner::createRunner($file, true);
    $runner = Tmdet30Runner::createRunner($file);
    // $runner->enableOverwrite = true;
    $runner->exec();
    // var_dump($runner->getCommandLine());die('END');
    echo "meld {$runner->oldTmdetFile} {$runner->newTmdetFile}" . PHP_EOL;
}
