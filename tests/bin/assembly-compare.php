<?php

use Unitmp\TmdetTest\Constants\FileSystem;
use Unitmp\TmdetTest\Services\Tmdet30Runner;

require_once 'vendor/autoload.php';

/*
 * This script compares PDBe assemblies with PDBTM 3.0 contained biomatricies.
 * It checks numbers of polymers in PDBe CIFs and PDBTM XMLs.
 */

const WORKING_DIR = 'issue_996';
if (!chdir(WORKING_DIR)) {
    fprintf(STDERR, "Couldn't change current directory to '%s'", WORKING_DIR);
    exit(1);
}
if (!file_exists(AssemblyCompare::XML_DIR)) {
    fprintf(STDERR, "Cache directory does not exist '%s'", AssemblyCompare::XML_DIR);
}

AssemblyCompare::main();

class AssemblyCompare {

    const PBDE_ASSEMBLY_URL = 'https://www.ebi.ac.uk/pdbe/static/entry/%s-assembly.xml';
    const XML_DIR = 'assembly-xmls'; // cache directory for xml files
    const HTTP_REQ_DELAY = 4; // 4 secs request delay

    public static function main() {

        // INPUT v1 - NOTE it results codeItems
        // $codes = array_slice(static::getPdbtmCodes(), 0, 40);
        //$codes = static::getPdbtmCodes();

        // INPUT v2 - NOTE it results only entry codes
        $diffResult = json_decode(file_get_contents('assembly-diff-v2.json'));
        $codes = array_merge($diffResult->does_not_match, $diffResult->match_but_not_first_assembly);
        // $codes = array_slice($codes, 0, 10);

        // differing assemblies
        $result = [
            'does_not_match' => [],
            'match_but_not_first_assembly' => [],
            'correctable_by_chain_delete' => [],
            'unexpected_failure' => [],
        ];
        foreach ($codes as $codeItem) {
            try {
                $code = $codeItem; //$codeItem['code'];
                $pdbtmXmlData = static::getPdbtmXmlData($code);
                $pdbtmChainCount = count($pdbtmXmlData->CHAIN);
                $pdbtmDeletedChainCount = static::getPdbtmDeletedChainCount($pdbtmXmlData);
                $assemblyData = static::getPdbEuropeAssemblyData($code);
                $pdbeChainCounts = static::getPdbEuropeChainCounts($assemblyData);

                $hasMatch = false;
                $firstMatches = false;
                foreach ($pdbeChainCounts as $assemblyItem) {
                    $correctableByChainDelete = $assemblyItem['chain_count'] == ($pdbtmChainCount + $pdbtmDeletedChainCount);
                    if ($assemblyItem['chain_count'] != $pdbtmChainCount || !$correctableByChainDelete) {
                        // skip this item
                        continue;
                    }

                    $hasMatch = true;
                    $asm = $assemblyItem['assembly_id'];
                    $detail = "";
                    if ($asm == "1") {
                        $firstMatches = true;
                    } else {
                        $detail = " - INFO: this is not assembly 1";
                    }
                    echo "$code: assembly$asm chain count is matching{$detail}\n";
                }

                if ($correctableByChainDelete) {
                    $result['correctable_by_chain_delete'][] = $code;
                    echo "$code: ERROR: NO MATCHING ASSEMBLY, MAYBE CORRECTABLE - expected chain count $pdbtmChainCount\n";
                    continue;
                }

                if (!$hasMatch) {
                    $result['does_not_match'][] = $code;
                    echo "$code: ERROR: THERE IS NO MATCHING ASSEMBLY - expected chain count $pdbtmChainCount\n";
                }

                if ($hasMatch && !$firstMatches) {
                    $result['match_but_not_first_assembly'][] = $code;
                    echo "$code: WARNING: DOES NOT MATCH WITH ASSEMBLY 1 - expected chain count $pdbtmChainCount\n";
                }

            } catch (Error|RuntimeException $e) {
                $code = $codeItem; //$codeItem['code'];
                $result['unexpected_failure'][] = $code;
                echo "ERROR: $code: {$e->getMessage()}" . PHP_EOL . $e->getTraceAsString() . PHP_EOL;
            }
        }
        echo "JSON RESULT:\n";
        echo "===================================\n";
        echo json_encode($result) . PHP_EOL;
    }

    public static function getPdbtmCodes(): array {
        $codes = Tmdet30Runner::getPdbtmCodes(FileSystem::PDBTM_ALL_XML_FILE);
        // return only TMP codes
        return array_filter($codes, function($value) { return $value['isTmp']; });
    }

    public static function getPdbEuropeAssemblyData($code): SimpleXMLElement {
        // iterate on pdbtm entries and their related PDBe assemblies
        $filePath = static::XML_DIR . "/$code-assembly.xml";
        if (file_exists($filePath)) {
            $url = $filePath;
        } else {
            $url = sprintf(static::PBDE_ASSEMBLY_URL, $code);
            sleep(static::HTTP_REQ_DELAY);
        }
        $content = file_get_contents($url);
        if (!file_exists($filePath)) {
            // cache the result, if it was donwloaded first time
            file_put_contents($filePath, $content);
        }
        return simplexml_load_string($content);
    }

    public static function getPdbEuropeChainCounts(SimpleXMLElement $assemblyData): array {

        $result = [];
        foreach($assemblyData->assembly as $assembly) {
            // echo "{$assembly['id']}; name: " . $assembly['name'] . PHP_EOL;
            $chainCount = 0;
            foreach ($assembly->entity as $entity) {
                if ($entity['type'] == 'polymer') {
                    // echo $entity['entity_id'] . ' ' . $entity['chain_ids'] . '; ch count: ' .  $entity['count'] . PHP_EOL;
                    $chainCount += intval($entity['count']);
                }
            }
            $result[] = [
                'assembly_id' => $assembly['id'],
                'chain_count' => $chainCount
            ];
        }
        return $result;
    }

    public static function getPdbtmXmlData(string $code): SimpleXMLElement {

        $subDir = '/' . $code[1] . $code[2] . '/';
        $fileName = "$code.xml";
        $fullPath = FileSystem::PDBTM30_DATA_DIR . "/database{$subDir}{$fileName}";
        return simplexml_load_file($fullPath);
    }

    public static function getPdbtmDeletedChainCount(SimpleXMLElement $xmlData): int {
        if ($xmlData->BIOMATRIX == null || $xmlData->BIOMATRIX->DELETE == null) {
            return 0;
        }
        return count($xmlData->BIOMATRIX->DELETE);
    }
}
