<?php

use Spatie\SimpleExcel\SimpleExcelWriter;
use Unitmp\TmdetTest\PdbtmDifferencesTest;
use Unitmp\TmdetTest\Services\DifferenceFlags;
use Unitmp\TmdetTest\Services\PdbtmComparator;
use Unitmp\TmdetTest\Services\Tmdet30Runner;

require_once 'vendor/autoload.php';

const COMPARISION_RESULT_FILE = './tmdet_differences.json';

// Pre-check write of output file
$outputFile = ($argc > 1) ? $argv[1] : COMPARISION_RESULT_FILE;
if (!file_put_contents($outputFile, "This file is writable by tmdiff. Waiting for JSON data...\n")) {
    fprintf(STDOUT, "Error: could not write '$outputFile'\n");
    exit(1);
}

TmDiff::main($outputFile);

class TmDiff {
    public static function main(string $outputFile) {

        $codes = PdbtmDifferencesTest::getEntPathsOfTmps();
        //$codes = array_slice($codes, 0, 1);
        $compResults = [];
        $byCategories = [];
        $fatalProcessingErrors = [];
        foreach ($codes as $file => $item) {
            echo "Processing data of $file" . PHP_EOL;
            $runner = Tmdet30Runner::createRunner($file);
            try {
                $runner->exec();

                echo "meld {$runner->oldTmdetFile} {$runner->newTmdetFile}" . PHP_EOL;
                $differences = PdbtmComparator::compareTmdetData($runner);

                $flags = $differences['flags'];
                // TODO
                // This condition is not strict enough yet.
                // It would be better to check presence/absence at chain level.
                if ($flags->onlyNewTmdetXmlHasDeletedChain) {
                    $runner = Tmdet30Runner::createRunner(entFile: $file, noForceDelete: true);
                    $runner->enableOverwrite = true;
                    $runner->noForceDelete = true;
                    echo "{$runner->pdbCode}: Running again with -no_forcedel" . PHP_EOL;
                    $runner->exec();
                }

            } catch (Exception $e) {
                $fatalProcessingErrors[] = $runner->pdbCode;
                fprintf(STDERR, "EXCEPTION: %s - %s\n", $runner->pdbCode, $e->getMessage());
                continue;
            }
            if (empty($differences['categories'])) {
                continue;
            }
            $result = [
                'code' => $runner->pdbCode,
                'details' => $differences
            ];
            $compResults[$runner->pdbCode] = $result;
            static::addCodeToCategories($byCategories, $result);
        }

        $byCategories['FATAL ERRORS DURING COMPARE PROCESS'] = $fatalProcessingErrors;
        ksort($byCategories);
        $byCategories = static::refactorCategories($byCategories);
        $outputData = [
            'byCategories' => $byCategories,
            'detailsByCodes' => $compResults
        ];

        echo "Writing result into $outputFile ... ";
        file_put_contents($outputFile,
            // json_encode($outputData, JSON_PRETTY_PRINT|JSON_UNESCAPED_SLASHES));
            json_encode($outputData, JSON_UNESCAPED_SLASHES));
        echo 'DONE.' . PHP_EOL;
        $outputFile = str_replace('.json', '.xlsx', $outputFile, $count);
        if ($count == 0) {
            $outputFile = "$outputFile.xlsx";
        }
        static::writeExcelOutput($outputData, "$outputFile");

    }

    public static function writeExcelOutput(array $compareResults, string $outputFile) {
        echo "Writing result in Excel format: $outputFile ... ";
        $writer = SimpleExcelWriter::create($outputFile);
        static::addOverviewSheet($writer, $compareResults);

        $categories = $compareResults['byCategories'];
        foreach ($categories as $index => $category) {
            $categoryName = $category['category'];
            $writer->nameCurrentSheet(substr($categoryName, 0, 31));
            $column = $categoryName;
            $writer->addHeader([ $column ]);
            foreach ($category['codes'] as $code) {
                $writer->addRow([ $column => $code ]);
            }
            if (count($categories) > $index + 1) {
                $writer->addNewSheetAndMakeItCurrent();
            }
        }
        $writer->close();
        echo 'DONE.' . PHP_EOL;
    }

    public static function addCodeToCategories(array &$categories, array &$compareResults): void {
        $code = $compareResults['code'];
        foreach ($compareResults['details']['categories'] as $category) {
            if (array_key_exists($category, $categories)) {
                $categories[$category][] = $code;
            } else {
                $categories[$category] = [ $code ];
            }
        }
    }

    public static function addOverviewSheet(SimpleExcelWriter $writer, array $results): void {

        $writer->nameCurrentSheet('Overview');

        // Construct header
        $reflection = new ReflectionClass(new DifferenceFlags());
        $columns = [ 'code' ];
        foreach ($reflection->getProperties() as $property) {
            $columns[] = $property->getName();
        }
        $writer->addHeader($columns);

        foreach ($results['detailsByCodes'] as $entry) {
            $flags = $entry['details']['flags'];
            $row = [ 'code' => $entry['code'] ];
            foreach ($flags as $key => $value) {
                $row[$key] = (int)$value;
            }
            $writer->addRow($row);
        }

        $writer->addNewSheetAndMakeItCurrent();
    }

    public static function refactorCategories(array &$categories): array {
        $refactored = [];
        foreach ($categories as $category => $items) {
            ksort($items);
            $refactored[] = [
                'category' => $category,
                'count' => count($items),
                'codes' => $items
            ];
        }
        return $refactored;
    }

}
