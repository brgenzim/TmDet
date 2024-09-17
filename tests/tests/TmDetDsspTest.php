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
    const MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT = 2;

    #[Group('dssp')]
    // #[DataProvider('cifPathsForAlignmentDebug')]
    #[DataProvider('cifPathsForSignificantDsspErrorTests')]
    // #[DataProvider('cifPathsWithErrorProvider')]
    // #[DataProvider('cifPathsForLenvenshteinFailuresDebug')]
    public function test_specific_dssp(string $cifPath) {
        // Pre-check
        $filter = new StructurePreFilter($cifPath);
        $exception = null;
        // if (!$filter->checkChainLengthsForDsspCalc()) {
        //     printf("Skipping test - not DSSP-eligible: %s\n", $cifPath);
        //     $this->markTestSkipped('There is a chain with unsupported length ' . $cifPath);
        // }
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

    //#[DataProvider('getCifPaths')]
    #[DataProvider('getCifPathsOfTmps')]
    public function test_whole_archive(string $cifPath) {

        // Pre-check
        $filter = new StructurePreFilter($cifPath);
        $exception = null;
        // if (!$filter->checkChainLengthsForDsspCalc()) {
        //     printf("Skipping test - not DSSP-eligible: %s\n", $cifPath);
        //     $this->markTestSkipped('There is a chain with unsupported length ' . $cifPath);
        // }

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

    #[DataProvider('getCifPathsOfTmpsWithErrors')]
    public function test_tmps_with_errors(string $cifPath) {

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
            $expectedOriginal = $expectedDssp;
            $expectedDssp = $this->promotifCompatible($expectedOriginal);
            $actualOriginal = $actualDssp = $actual[$key];
            $actualDssp = $this->promotifCompatible($actualDssp);
            if (strlen($expectedDssp) != strlen($actualDssp)) {
                break;
            }
            // If old tmdet does not recognize DSSP string
            if (preg_match('/^-+$/', $expectedDssp, $matches)) {
                continue;
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

        if (!isset($expectedDssp)) {
            $this->assertTrue(true);
            return;
        }

        $this->assertEquals(strlen($expectedDssp), strlen($actualDssp),
            'DSSP compare failed at chain ' . $key . ' - dssp string lengths differ'
            . "\nexpected: '$expectedDssp'"
            . "\nactual:   '$actualDssp'"
            . "\nOriginal DSSP string before 'promotif' conversion:"
            . "\nexpected: '$expectedOriginal'"
            . "\nactual:   '$actualOriginal'");

        $this->assertTrue($ok, 'DSSP compare failed at chain ' . $key
            . "\nexpected: '$expectedDssp'"
            . "\nactual:   '$actualDssp'"
            . "\nexpected: '$expectedOriginal'"
            . "\nactual:   '$actualOriginal'"
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
        return [ [ 'not-existing.cif.gz' ] ];

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
        return array_slice($files, 41000, 11000);
        //return $files;
    }

    public static function getCifPathsOfTmps(): array {
        $cifPaths = [];
        foreach (explode("\n", trim(file_get_contents('../data/tm-codes.txt'))) as $code) {
            $fileName = "$code.cif.gz";
            $subDir = '/' . $code[1] . $code[2] . '/';
            $cifPaths[$fileName] = [ static::PDB_ZFS_DIR . $subDir . $fileName  ];
        }
        return $cifPaths;
    }

    public static function getCifPathsOfTmpsWithErrors(): array {
        return [
            '1bcc.cif.gz' => [ static::PDB_ZFS_DIR . '/bc/1bcc.cif.gz' ],
            '7a4p.cif.gz' => [ static::PDB_ZFS_DIR . '/a4/7a4p.cif.gz' ],
            '7t6b.cif.gz' => [ static::PDB_ZFS_DIR . '/t6/7t6b.cif.gz' ],
            '8e4o.cif.gz' => [ static::PDB_ZFS_DIR . '/e4/8e4o.cif.gz' ],
            '3bcc.cif.gz' => [ static::PDB_ZFS_DIR . '/bc/3bcc.cif.gz' ],
            '4msw.cif.gz' => [ static::PDB_ZFS_DIR . '/ms/4msw.cif.gz' ],
            '6zxs.cif.gz' => [ static::PDB_ZFS_DIR . '/zx/6zxs.cif.gz' ],
            '7a5v.cif.gz' => [ static::PDB_ZFS_DIR . '/a5/7a5v.cif.gz' ],
            '2axt.cif.gz' => [ static::PDB_ZFS_DIR . '/ax/2axt.cif.gz' ],
            '2bcc.cif.gz' => [ static::PDB_ZFS_DIR . '/bc/2bcc.cif.gz' ],
        ];
    }

    /**
     * For issue https://redmine.enzim.ttk.hu/issues/895
     */
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
            '8e74.cif.gz' => [ static::PDB_ZFS_DIR . '/e7/8e74.cif.gz' ],
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

    public static function cifPathsForLenvenshteinFailuresDebug(): array {
        return [
            '5ei1.cif.gz' => [ static::PDB_ZFS_DIR . '/ei/5ei1.cif.gz' ],
            '5eit.cif.gz' => [ static::PDB_ZFS_DIR . '/ei/5eit.cif.gz' ],
            '8el5.cif.gz' => [ static::PDB_ZFS_DIR . '/el/8el5.cif.gz' ],
            '5eqq.cif.gz' => [ static::PDB_ZFS_DIR . '/eq/5eqq.cif.gz' ],
            '8ev2.cif.gz' => [ static::PDB_ZFS_DIR . '/ev/8ev2.cif.gz' ],
            '4exw.cif.gz' => [ static::PDB_ZFS_DIR . '/ex/4exw.cif.gz' ],
            '8f0l.cif.gz' => [ static::PDB_ZFS_DIR . '/f0/8f0l.cif.gz' ],
            '5f1w.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/5f1w.cif.gz' ],
            '6f1x.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/6f1x.cif.gz' ],
            '8f12.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/8f12.cif.gz' ],
            '8f14.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/8f14.cif.gz' ],
            '8f17.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/8f17.cif.gz' ],
            '7f22.cif.gz' => [ static::PDB_ZFS_DIR . '/f2/7f22.cif.gz' ],
            '6f9w.cif.gz' => [ static::PDB_ZFS_DIR . '/f9/6f9w.cif.gz' ],
            '5ef5.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/5ef5.cif.gz' ],
            '4egy.cif.gz' => [ static::PDB_ZFS_DIR . '/eg/4egy.cif.gz' ],
            '3eh4.cif.gz' => [ static::PDB_ZFS_DIR . '/eh/3eh4.cif.gz' ],
            '8el4.cif.gz' => [ static::PDB_ZFS_DIR . '/el/8el4.cif.gz' ],
            '3et1.cif.gz' => [ static::PDB_ZFS_DIR . '/et/3et1.cif.gz' ],
            '7f20.cif.gz' => [ static::PDB_ZFS_DIR . '/f2/7f20.cif.gz' ],
            '3f43.cif.gz' => [ static::PDB_ZFS_DIR . '/f4/3f43.cif.gz' ],
            '1ee7.cif.gz' => [ static::PDB_ZFS_DIR . '/ee/1ee7.cif.gz' ],
            '2efg.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/2efg.cif.gz' ],
            '3eh3.cif.gz' => [ static::PDB_ZFS_DIR . '/eh/3eh3.cif.gz' ],
            '5ehj.cif.gz' => [ static::PDB_ZFS_DIR . '/eh/5ehj.cif.gz' ],
            '8el3.cif.gz' => [ static::PDB_ZFS_DIR . '/el/8el3.cif.gz' ],
            '8f0z.cif.gz' => [ static::PDB_ZFS_DIR . '/f0/8f0z.cif.gz' ],
            '6f1y.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/6f1y.cif.gz' ],
            '8f16.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/8f16.cif.gz' ],
            '3f5h.cif.gz' => [ static::PDB_ZFS_DIR . '/f5/3f5h.cif.gz' ],
            '6f5j.cif.gz' => [ static::PDB_ZFS_DIR . '/f5/6f5j.cif.gz' ],
            '7ec3.cif.gz' => [ static::PDB_ZFS_DIR . '/ec/7ec3.cif.gz' ],
            '6een.cif.gz' => [ static::PDB_ZFS_DIR . '/ee/6een.cif.gz' ],
            '5egv.cif.gz' => [ static::PDB_ZFS_DIR . '/eg/5egv.cif.gz' ],
            '3eh5.cif.gz' => [ static::PDB_ZFS_DIR . '/eh/3eh5.cif.gz' ],
            '8el6.cif.gz' => [ static::PDB_ZFS_DIR . '/el/8el6.cif.gz' ],
            '8eq6.cif.gz' => [ static::PDB_ZFS_DIR . '/eq/8eq6.cif.gz' ],
            '8ev1.cif.gz' => [ static::PDB_ZFS_DIR . '/ev/8ev1.cif.gz' ],
            '6exv.cif.gz' => [ static::PDB_ZFS_DIR . '/ex/6exv.cif.gz' ],
            '8f10.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/8f10.cif.gz' ],
            '8f13.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/8f13.cif.gz' ],
            '8f15.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/8f15.cif.gz' ],
            '7f21.cif.gz' => [ static::PDB_ZFS_DIR . '/f2/7f21.cif.gz' ],
            '7f4g.cif.gz' => [ static::PDB_ZFS_DIR . '/f4/7f4g.cif.gz' ],
            '6f5u.cif.gz' => [ static::PDB_ZFS_DIR . '/f5/6f5u.cif.gz' ],
            '6f7w.cif.gz' => [ static::PDB_ZFS_DIR . '/f7/6f7w.cif.gz' ],
            '7f7g.cif.gz' => [ static::PDB_ZFS_DIR . '/f7/7f7g.cif.gz' ],
        ];
    }

    public static function cifPathsForAlignmentDebug(): array {
        return [
            '4em2.cif.gz' => [ static::PDB_ZFS_DIR . '/em/4em2.cif.gz' ],
            '1afo.cif.gz' => [ static::PDB_ZFS_DIR . '/af/1afo.cif.gz' ],
            '5eit.cif.gz' => [ static::PDB_ZFS_DIR . '/ei/5eit.cif.gz' ],
            '8f14.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/8f14.cif.gz' ],
            '3et1.cif.gz' => [ static::PDB_ZFS_DIR . '/et/3et1.cif.gz' ],
            '7f7g.cif.gz' => [ static::PDB_ZFS_DIR . '/f7/7f7g.cif.gz' ],
        ];
    }

    public static function cifPathsForDsspLengthErrorDebug(): array {

        return [
            '6ecs.cif.gz' => [ static::PDB_ZFS_DIR . '/ec/6ecs.cif.gz' ],
            '7ega.cif.gz' => [ static::PDB_ZFS_DIR . '/eg/7ega.cif.gz' ],
            '8eih.cif.gz' => [ static::PDB_ZFS_DIR . '/ei/8eih.cif.gz' ],
            '6ejl.cif.gz' => [ static::PDB_ZFS_DIR . '/ej/6ejl.cif.gz' ],
            '6em3.cif.gz' => [ static::PDB_ZFS_DIR . '/em/6em3.cif.gz' ],
            '7etm.cif.gz' => [ static::PDB_ZFS_DIR . '/et/7etm.cif.gz' ],
            '8etc.cif.gz' => [ static::PDB_ZFS_DIR . '/et/8etc.cif.gz' ],
            '4ewc.cif.gz' => [ static::PDB_ZFS_DIR . '/ew/4ewc.cif.gz' ],
            '6exn.cif.gz' => [ static::PDB_ZFS_DIR . '/ex/6exn.cif.gz' ],
            '8f6n.cif.gz' => [ static::PDB_ZFS_DIR . '/f6/8f6n.cif.gz' ],
            '5fbt.cif.gz' => [ static::PDB_ZFS_DIR . '/fb/5fbt.cif.gz' ],
            '7edx.cif.gz' => [ static::PDB_ZFS_DIR . '/ed/7edx.cif.gz' ],
            '2efw.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/2efw.cif.gz' ],
            '2efz.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/2efz.cif.gz' ],
            '3ef1.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/3ef1.cif.gz' ],
            '7efd.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/7efd.cif.gz' ],
            '7efh.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/7efh.cif.gz' ],
            '7efm.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/7efm.cif.gz' ],
            '5eg2.cif.gz' => [ static::PDB_ZFS_DIR . '/eg/5eg2.cif.gz' ],
            '7eg9.cif.gz' => [ static::PDB_ZFS_DIR . '/eg/7eg9.cif.gz' ],
            '4em0.cif.gz' => [ static::PDB_ZFS_DIR . '/em/4em0.cif.gz' ],
            '4em2.cif.gz' => [ static::PDB_ZFS_DIR . '/em/4em2.cif.gz' ],
            '6em5.cif.gz' => [ static::PDB_ZFS_DIR . '/em/6em5.cif.gz' ],
            '7eo4.cif.gz' => [ static::PDB_ZFS_DIR . '/eo/7eo4.cif.gz' ],
            '8epl.cif.gz' => [ static::PDB_ZFS_DIR . '/ep/8epl.cif.gz' ],
            '8ev3.cif.gz' => [ static::PDB_ZFS_DIR . '/ev/8ev3.cif.gz' ],
            '8ezr.cif.gz' => [ static::PDB_ZFS_DIR . '/ez/8ezr.cif.gz' ],
            '7efg.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/7efg.cif.gz' ],
            '7efi.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/7efi.cif.gz' ],
            '7efn.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/7efn.cif.gz' ],
            '7eg8.cif.gz' => [ static::PDB_ZFS_DIR . '/eg/7eg8.cif.gz' ],
            '8eii.cif.gz' => [ static::PDB_ZFS_DIR . '/ei/8eii.cif.gz' ],
            '7ept.cif.gz' => [ static::PDB_ZFS_DIR . '/ep/7ept.cif.gz' ],
            '8eti.cif.gz' => [ static::PDB_ZFS_DIR . '/et/8eti.cif.gz' ],
            '8eup.cif.gz' => [ static::PDB_ZFS_DIR . '/eu/8eup.cif.gz' ],
            '6f09.cif.gz' => [ static::PDB_ZFS_DIR . '/f0/6f09.cif.gz' ],
            '6f2r.cif.gz' => [ static::PDB_ZFS_DIR . '/f2/6f2r.cif.gz' ],
            '3f5b.cif.gz' => [ static::PDB_ZFS_DIR . '/f5/3f5b.cif.gz' ],
            '2f69.cif.gz' => [ static::PDB_ZFS_DIR . '/f6/2f69.cif.gz' ],
            '8edo.cif.gz' => [ static::PDB_ZFS_DIR . '/ed/8edo.cif.gz' ],
            '3ee0.cif.gz' => [ static::PDB_ZFS_DIR . '/ee/3ee0.cif.gz' ],
            '7ef3.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/7ef3.cif.gz' ],
            '7ef9.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/7ef9.cif.gz' ],
            '7efl.cif.gz' => [ static::PDB_ZFS_DIR . '/ef/7efl.cif.gz' ],
            '7eg7.cif.gz' => [ static::PDB_ZFS_DIR . '/eg/7eg7.cif.gz' ],
            '6em4.cif.gz' => [ static::PDB_ZFS_DIR . '/em/6em4.cif.gz' ],
            '8euy.cif.gz' => [ static::PDB_ZFS_DIR . '/eu/8euy.cif.gz' ],
            '8f2k.cif.gz' => [ static::PDB_ZFS_DIR . '/f2/8f2k.cif.gz' ],
            '4f88.cif.gz' => [ static::PDB_ZFS_DIR . '/f8/4f88.cif.gz' ],
        ];
    }

    public static function cifPathsForSignificantDsspErrorTests(): array {
        return [
            '5eqq.cif.gz' => [ static::PDB_ZFS_DIR . '/eq/5eqq.cif.gz' ],
            '3f43.cif.gz' => [ static::PDB_ZFS_DIR . '/f4/3f43.cif.gz' ],
            '1ee7.cif.gz' => [ static::PDB_ZFS_DIR . '/ee/1ee7.cif.gz' ],
            '4exw.cif.gz' => [ static::PDB_ZFS_DIR . '/ex/4exw.cif.gz' ],
            '6f1x.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/6f1x.cif.gz' ],
            '4egy.cif.gz' => [ static::PDB_ZFS_DIR . '/eg/4egy.cif.gz' ],
            '3ee0.cif.gz' => [ static::PDB_ZFS_DIR . '/ee/3ee0.cif.gz' ],
            '3f5b.cif.gz' => [ static::PDB_ZFS_DIR . '/f5/3f5b.cif.gz' ],
            '6f5j.cif.gz' => [ static::PDB_ZFS_DIR . '/f5/6f5j.cif.gz' ],
            '6f7w.cif.gz' => [ static::PDB_ZFS_DIR . '/f7/6f7w.cif.gz' ],

            '4egy.cif.gz' => [ static::PDB_ZFS_DIR . '/eg/4egy.cif.gz' ],
            '7e99.cif.gz' => [ static::PDB_ZFS_DIR . '/e9/7e99.cif.gz' ],

            '5f1w.cif.gz' => [ static::PDB_ZFS_DIR . '/f1/5f1w.cif.gz' ],
        ];
    }
}
