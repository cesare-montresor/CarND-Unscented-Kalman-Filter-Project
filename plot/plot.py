import matplotlib.pyplot as plt
import csv
from scipy.stats import chi2

radar = {
    'name': 'Radar',
    'df': 3,
    'path': '../ide_profiles/xcode/Debug/nis_radar.csv'
}

lidar = {
    'name': 'Lidar',
    'df': 2,
    'path': '../ide_profiles/xcode/Debug/nis_lidar.csv'
}

sensors = [radar, lidar]

for sensor in sensors:
    name = sensor['name']
    path = sensor['path']
    df = sensor['df']

    lines = []
    with open(path, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        lines = list(csvreader)

    chi2_050 = chi2.isf(0.05, df)

    plt.axhline(chi2_050, color='r')
    plt.plot(lines)
    plt.title(name)
    plt.ylabel('NIS')
    plt.xlabel('Measurements')
    plt.legend(['X.050', 'NIS'], loc='upper right')
    plt.show()
    plt.close('all')
