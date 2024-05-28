<?php

namespace Unitmp\TmdetTest;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\Attributes\Group;
use PHPUnit\Framework\TestCase;
use Unitmp\TmdetTest\Constants\FileSystem;
use Unitmp\TmdetTest\Services\PdbtmDsspRunner;
use Unitmp\TmdetTest\Services\StructurePreFilter;
use Unitmp\TmdetTest\Services\TmDetDsspRunner;

class TmDetDsspTest extends TestCase {

    const PDB_ZFS_DIR = FileSystem::PDB_ZFS_DIR;
    const PDB_ENT_ZFS_DIR = FileSystem::PDB_ENT_ZFS_DIR;
    const MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT = 15;

    #[Group('dssp')]
    public function test_pdbtm_dssp_integrity_with_7rh7_chain_w() {
        // Arrange
        $pdbtmRunner = PdbtmDsspRunner::createRunner(static::PDB_ENT_ZFS_DIR . '/rh/pdb7rh7.ent.gz');
        $tmdetRunner = TmDetDsspRunner::createRunner(static::PDB_ZFS_DIR . '/rh/7rh7.cif.gz');
        // Act
        $isSuccess = $pdbtmRunner->exec();
        $isTmDetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($isSuccess);
        $this->assertTrue($isTmDetSuccess);
        $this->assertEquals($pdbtmRunner->dssps['W'], $tmdetRunner->dssps['W']);
        $this->assertEquals($pdbtmRunner->dssps, $tmdetRunner->dssps);
    }

    #[Group('dssp')]
    #[DataProvider('cifPathsWithErrorProvider')]
    public function test_specific_dssp(string $cifPath) {
        // Pre-check
        $filter = new StructurePreFilter($cifPath);
        $exception = null;
        if (!$filter->tryCheckEntryFiles($exception)) {
            printf("Skipping test: %s\n", $exception->getMessage());
            $this->markTestSkipped('Unsupported CIF-ENT pair ' . $cifPath);
        }

        // Arrange
        $tmdetRunner = TmDetDsspRunner::createRunner($cifPath);
        $dsspRunner = PdbtmDsspRunner::createRunner($filter->entFile);
        // Act
        $dsspSuccess = $dsspRunner->exec();
        $tmdetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($dsspSuccess);
        $this->assertTrue($tmdetSuccess);
        $this->assertDsspArraysEqual($dsspRunner->dssps, $tmdetRunner->dssps);
    }

    #[DataProvider('getCifPaths')]
    public function test_whole_archive(string $cifPath) {

        // Pre-check
        $filter = new StructurePreFilter($cifPath);
        $exception = null;
        if (!$filter->tryCheckEntryFiles($exception)) {
            printf("Skipping test: %s\n", $exception->getMessage());
            $this->markTestSkipped('Unsupported CIF-ENT pair ' . $cifPath);
        }

        // Arrange
        $tmdetRunner = TmDetDsspRunner::createRunner($cifPath);
        $dsspRunner = PdbtmDsspRunner::createRunner($filter->entFile);
        // Act
        $dsspSuccess = $dsspRunner->exec();
        $tmdetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertTrue($dsspSuccess);
        $this->assertDsspArraysEqual($dsspRunner->dssps, $tmdetRunner->dssps);
    }

    //
    // Util functions
    //

    public function assertDsspArraysEqual(array $expected, array $actual) {
        $this->assertEquals(array_keys($expected), array_keys($actual));
        $ok = true;
        foreach ($expected as $key => $expectedDssp) {
            $expectedDssp = $this->promotifCompatible($expectedDssp);
            $actualDssp = $actual[$key];
            if (strlen($expectedDssp) != strlen($actualDssp)) {
                break;
            }
            // max 2% difference allowed
            $maxDistance = round(strlen($expectedDssp) / 100 * static::MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT, 2);
            $actualDistance = levenshtein($expectedDssp, $actualDssp);
            $actualPercent = round($actualDistance / strlen($actualDssp) * 100, 2);
            if ($actualDistance > $maxDistance) {
                $ok = false;
                break;
            }
        }
        $this->assertEquals(strlen($expectedDssp), strlen($actualDssp),
            'DSSP compare failed at chain ' . $key
            . "\nexpected: '$expectedDssp'"
            . "\nactual:   '$actualDssp'");

        $this->assertTrue($ok, 'DSSP compare failed at chain ' . $key
            . "\nexpected: '$expectedDssp'"
            . "\nactual:   '$actualDssp'"
            . "\nmax expected levenshtein distance (%): "
            .  static::MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT
            . ", actual distance (%): $actualPercent"
            . "\nmax expected levenshtein distance: $maxDistance, actual distance: $actualDistance");
    }

    public function promotifCompatible(string $dssp): string {
        // reduce to PROMOTIF-compatible version (MAXIT: process_entry)
        $dssp = preg_replace('/[GIP]/', 'H', $dssp);

        // return preg_replace('/[SB]/', '-', $dssp);
        // in most cases T is set
        return preg_replace('/[STB]/', '-', $dssp);
    }

    /**
     * Data provider for test_whole_archive
     */
    public static function getCifPaths(): array {

        // TODO: remove later
        //return [ [ 'not-existing.cif.gz' ] ];

        printf("Scanning PDB directories...\n");
        $dirs = scandir(static::PDB_ZFS_DIR);
        if (!$dirs || count($dirs) == 2) {
            return [];
        }

        $files = [];
        foreach (array_slice($dirs, 2) as $directory) {
            $directory = static::PDB_ZFS_DIR . "/$directory";
            $items = scandir($directory);
            $items = preg_grep('/^.+cif.gz$/', $items);
            foreach($items as $item) {
                $files[ $item ] = [ "$directory/$item" ];
            }
        }
        // TODO: remove slice later
        //return array_slice($files, 100, 25);
        //return array_slice($files, 100, 1000);
        return array_slice($files, 40000, 1000);
        //return $files;
    }

    public static function cifPathsWithErrorProvider(): array {
        return [
            '5e8d.cif.gz' => [ static::PDB_ZFS_DIR . '/e8/5e8d.cif.gz' ],
            '3eap.cif.gz' => [ static::PDB_ZFS_DIR . '/ea/3eap.cif.gz' ],
            '8eay.cif.gz' => [ static::PDB_ZFS_DIR . '/ea/8eay.cif.gz' ],
            '6ebd.cif.gz' => [ static::PDB_ZFS_DIR . '/eb/6ebd.cif.gz' ],
            '6ebq.cif.gz' => [ static::PDB_ZFS_DIR . '/eb/6ebq.cif.gz' ],
            '5e9d.cif.gz' => [ static::PDB_ZFS_DIR . '/e9/5e9d.cif.gz' ],
            '6eb4.cif.gz' => [ static::PDB_ZFS_DIR . '/eb/6eb4.cif.gz' ],
            '6ebg.cif.gz' => [ static::PDB_ZFS_DIR . '/eb/6ebg.cif.gz' ],
            '6e8r.cif.gz' => [ static::PDB_ZFS_DIR . '/e8/6e8r.cif.gz' ],
            '6e9q.cif.gz' => [ static::PDB_ZFS_DIR . '/e9/6e9q.cif.gz' ],
            '7e99.cif.gz' => [ static::PDB_ZFS_DIR . '/e9/7e99.cif.gz' ],
            '4e8o.cif.gz' => [ static::PDB_ZFS_DIR . '/e8/4e8o.cif.gz' ],
            '4e98.cif.gz' => [ static::PDB_ZFS_DIR . '/e9/4e98.cif.gz' ],
            '6ebc.cif.gz' => [ static::PDB_ZFS_DIR . '/eb/6ebc.cif.gz' ],
           
        ];
    }

    public static function regressionCifPathProvider(): array {
        return [
            '7rh7.cif.gz' => [ static::PDB_ZFS_DIR . '/rh/7rh7.cif.gz' ],
            '6a1t.cif.gz' => [ static::PDB_ZFS_DIR . '/a1/6a1t.cif.gz' ],
            '7eb2.cif.gz' => [ static::PDB_ZFS_DIR . '/eb/7eb2.cif.gz' ],
            '1e8x.cif.gz' => [ static::PDB_ZFS_DIR . '/e8/1e8x.cif.gz' ],
        ];
    }

}
