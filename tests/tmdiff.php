<?php

use Spatie\SimpleExcel\SimpleExcelWriter;
use Unitmp\TmdetTest\PdbtmDifferencesTest;
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
        // $codes = array_slice($codes, 0, 10);
        $compResults = [];
        $byCategories = [];
        $fatalProcessingErrors = [];
        foreach ($codes as $file => $item) {
            echo "Processing data of $file" . PHP_EOL;
            $runner = Tmdet30Runner::createRunner($file);
            // $runner->enableOverwrite = true;
            try {
                $runner->exec();
            } catch (Exception $e) {
                $fatalProcessingErrors[] = $runner->pdbCode;
                fprintf(STDERR, "EXCEPTION: %s - %s\n", $runner->pdbCode, $e->getMessage());
                continue;
            }
            echo "meld {$runner->oldTmdetFile} {$runner->newTmdetFile}" . PHP_EOL;
            $differences = PdbtmComparator::compareTmdetData($runner);
            if (empty($differences['categories'])) {
                continue;
            }
            $result = [
                'code' => $runner->pdbCode,
                'details' => $differences
            ];
            $compResults[$runner->pdbCode] = $result;
            static::addCodeToCategory($byCategories, $result);
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
        $categories = $compareResults['byCategories'];
        foreach ($categories as $index => $category) {
            $categoryName = $category['category'];
            $writer->nameCurrentSheet(substr($categoryName, 0, 31));
            $column = 'PDB Codes (' . $categoryName . ')';
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

    public static function addCodeToCategory(array &$categories, array &$compareResults): void {
        $code = $compareResults['code'];
        foreach ($compareResults['details']['categories'] as $category) {
            if (array_key_exists($category, $categories)) {
                $categories[$category][] = $code;
            } else {
                $categories[$category] = [ $code ];
            }
        }
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
