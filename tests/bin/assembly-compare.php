<?php

use Unitmp\TmdetTest\Constants\FileSystem;
use Unitmp\TmdetTest\Services\Tmdet30Runner;

require_once 'vendor/autoload.php';

/*
 * This script compares PDBe assemblies with PDBTM 3.0 contained biomatricies.
 * It checks numbers of polymers in PDBe CIFs and PDBTM XMLs.
 */

AssemblyCompare::main();

class AssemblyCompare {

    const PBDE_ASSEMBLY_URL = 'https://www.ebi.ac.uk/pdbe/static/entry/%s-assembly.xml';
    const HTTP_REQ_DELAY = 4; // 4 secs request delay

    public static function main() {

        // $codes = array_slice(static::getPdbtmCodes(), 0, 40);
        $codes = static::getPdbtmCodes();

        // differing assemblies
        $result = [
            'does_not_match' => [],
            'match_but_not_first_assembly' => [],
            'unexpected_failure' => [],
        ];
        foreach ($codes as $codeItem) {
            try {
                $code = $codeItem['code'];
                $pdbtmChainCount = static::getPdbtmChainCountOfBioMx($code);
                $pdbeChainCounts = static::getPdbEuropeChainCounts($code);

                $hasMatch = false;
                $firstMatches = false;
                foreach ($pdbeChainCounts as $assemblyItem) {
                    if ($assemblyItem['chain_count'] != $pdbtmChainCount) {
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

                if (!$hasMatch) {
                    $result['does_not_match'][] = $code;
                    echo "$code: ERROR: THERE IS NO MATCHING ASSEMBLY - expected chain count $pdbtmChainCount\n";
                }

                if ($hasMatch && !$firstMatches) {
                    $result['match_but_not_first_assembly'][] = $code;
                    echo "$code: WARNING: DOES NOT MATCH WITH ASSEMBLY 1 - expected chain count $pdbtmChainCount\n";
                }

                sleep(static::HTTP_REQ_DELAY);
            } catch (Error|RuntimeException $e) {
                $code = $codeItem['code'];
                $result['unexpected_failure'][] = $code;
                echo "ERROR: $code: " . $e->getTraceAsString() . PHP_EOL;
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

    public static function getPdbEuropeChainCounts($code): array {
        // iterate on pdbtm entries and their related PDBe assemblies
        $url = sprintf(static::PBDE_ASSEMBLY_URL, $code);

        $result = [];
        foreach(simplexml_load_file($url)->assembly as $assembly) {
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

    public static function getPdbtmChainCountOfBioMx(string $code): int {

        $subDir = '/' . $code[1] . $code[2] . '/';
        $fileName = "$code.xml";
        $fullPath = FileSystem::PDBTM30_DATA_DIR . "/database{$subDir}{$fileName}";
        $xmlObj = simplexml_load_file($fullPath);

        return count($xmlObj->CHAIN);
    }

}
