from galaxy import galaxy


def test_galaxy_init():
    MW = galaxy.Galaxy('MW', snap=0)
    assert MW.filename == 'MW_000.txt'

MW = galaxy.Galaxy('MW', snap=0)

def test_MW_data():
    assert MW.data.shape == (135000, )

def test_filter_by_type():
    type = 2 # just disk particles
    disk_particles = MW.filter_by_type(2)
    assert disk_particles.shape == (75000,)