<?php

namespace Unitmp\TmdetTest;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\Attributes\Group;
use PHPUnit\Framework\TestCase;
use Unitmp\TmdetTest\Constants\FileSystem;
use Unitmp\TmdetTest\Services\Tmdet30Runner;

use function PHPUnit\Framework\assertEquals;

class PdbtmDifferencesTest extends TestCase {

    const PDB_ENT_DIR = FileSystem::PDB_ENT_ZFS_DIR;
    const PDBTM_ALL_XML_FILE = FileSystem::PDBTM_ALL_XML_FILE;
    const MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT = 40;
    const MAX_ALLOWED_RESIDUE_DIFFERENCE = 4;

    #[Group('tmdet')]
    #[DataProvider('getEntPathsOfTmps')]
    public function test_pdbtm_manual_updates(string $entPath) {
        // Arrange
        $tmdetRunner = Tmdet30Runner::createRunner($entPath);
        if ($tmdetRunner->oldData['tmType'] !== 'Tm_Alpha') {
            printf("Skipping test (not Tm_Alpha): %s\n", $tmdetRunner->oldTmdetFile);
            $this->markTestSkipped('Not Tm_Alpha ' . $entPath);
        }

        // Act
        $tmdetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertTmpValueEquals($tmdetRunner);
        //$this->assertRegionArraysEqual($tmdetRunner->oldData, $tmdetRunner->newData);
    }

    //
    // Util functions
    //

    /**
     * Assert TM protein values are equal or their difference is small enough.
     */
    public function assertTmpValueEquals(Tmdet30Runner $runner) {
        $oldValue = $runner->oldData['isTmp'];
        $newValue = $runner->newData['isTmp'];
        $this->assertEquals($oldValue, $newValue, '"isTmp" values differ');

        $oldValue = $runner->oldData['chains'];
        $newValue = $runner->newData['chains'];
        $this->assertEquals(array_keys($oldValue), array_keys($newValue),
            'Different chain lists in TMDET results');

        $helixFilter = function ($regions) {
            return array_filter($regions, function ($value) {
                return $value === 'H';
            } );
        };
        foreach ($oldValue as $chain => $chainData) {
            assertEquals($chainData['num_tm'], $newValue[$chain]['num_tm'],
                'Number of TM regions differs');
            $newRegions = $helixFilter($newValue[$chain]['regions']);
            foreach ($helixFilter($chainData['regions']) as $index => $oldRegion) {
                $newRegion = $newRegions[$index];
                if ($oldRegion['type'] == 'H') {
                    $this->assertEquals($oldRegion['type'], $newRegion['type'],
                        'Region type of chain "' . $chain . '"[new seq_beg: ' . $newRegion['seq_beg'] . '] differs');
                }
                $this->assertEqualsWithDelta($oldRegion['seq_beg'], $newRegion['seq_beg'], static::MAX_ALLOWED_RESIDUE_DIFFERENCE,
                    'Region start in chain "' . $chain . '"[new seq_beg: ' . $newRegion['seq_beg'] . '] differs');
                $this->assertEqualsWithDelta($oldRegion['seq_end'], $newRegion['seq_end'], static::MAX_ALLOWED_RESIDUE_DIFFERENCE,
                    'Region end in chain "' . $chain . '"[new seq_end: ' . $newRegion['seq_end'] . '] differs');
            }
        }
    }

    public function assertRegionArraysEqual(array $expected, array $actual) {
        $this->assertEquals(array_keys($expected), array_keys($actual));
        $ok = true;
        foreach ($expected as $key => $expectedRegionString) {
            $actualRegionString = $actual[$key];
            if (strlen($expectedRegionString) != strlen($actualRegionString)) {
                break;
            }
            // max 2% difference allowed
            $maxDistance = round(strlen($expectedRegionString) / 100 * static::MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT, 2);
            $actualDistance = levenshtein($expectedRegionString, $actualRegionString);
            $actualPercent = round($actualDistance / strlen($actualRegionString) * 100, 2);
            if ($actualDistance > $maxDistance) {
                $ok = false;
                break;
            }
        }

        if (!isset($expectedRegionString)) {
            $this->assertTrue(true);
            return;
        }

        $this->assertEquals(strlen($expectedRegionString), strlen($actualRegionString),
            'Region compare failed at chain ' . $key . ' - region string lengths differ'
            . "\nexpected: '$expectedRegionString'"
            . "\nactual:   '$actualRegionString'");

        $this->assertTrue($ok, 'Region compare failed at chain ' . $key
            . "\nexpected: '$expectedRegionString'"
            . "\nactual:   '$actualRegionString'"
            . "\nmax expected levenshtein distance (%): "
            .  static::MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT
            . ", actual distance (%): $actualPercent"
            . "\nmax expected levenshtein distance: $maxDistance, actual distance: $actualDistance");
    }

    public static function getEntPathsOfTmps(): array {
        $codes = Tmdet30Runner::getPdbtmCodes(static::PDBTM_ALL_XML_FILE);

        $paths = [];
        foreach ($codes as $code => $item) {
            $subDir = '/' . $code[1] . $code[2] . '/';
            $fileName = "pdb$code.ent.gz";
            $fullPath = static::PDB_ENT_DIR . $subDir . $fileName;
            $paths[ $fileName ] = [ $fullPath ]; //  , $item['isTmp']
        }
        ksort($paths);
        // return array_slice($paths, 0, 100);
        return $paths;
    }

}
