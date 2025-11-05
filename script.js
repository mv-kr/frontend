// Ensure the DOM is fully loaded before accessing elements
let data = [];
let currentChart = null; 
let currentCharts = []; 
let currentSettings = null;
let dataset_id = null;
let uniqueNames = new Map();
// Define the settings to display and their descriptive labels
const settingsToDisplay = {
  'general': ['identifier_str', 'threads_no', 'chain_repeats', 'seed'],
  'MCMC_gibbs_sMMALA': ['iterations', 'iter_adapt'],
  'shMCMC_gibbs_sMMALA': ['iterations','adapt_runs']
};

const groupLabels = {
  'general': 'General Settings',
  'MCMC_gibbs_sMMALA': 'MCMC Settings',
  'shMCMC_gibbs_sMMALA': 'short MCMC'
};

const settingLabels = {
  'general': {
    'identifier_str': 'Job Identifier',
    'threads_no': 'Number of Threads',
    'chain_repeats': 'Chain Repeats',
    'seed': 'Seed'
  },
  'MCMC_gibbs_sMMALA': {
    'iterations': 'Number of Iterations',
    'iter_adapt': 'Adaptation Period Iterations'
  },
  'shMCMC_gibbs_sMMALA': {
    'iterations': 'Number of Iterations',
    'adapt_runs': 'Number of Adaptation Runs'
  }
};



const findOutputFolderBtn = document.getElementById('find-output-folder');
const resultsOutputInput = document.getElementById('results-output');

findOutputFolderBtn.addEventListener('click', () => {
    window.electron.openFolderDialog()
      .then(result => {
        if (!result.canceled) {
          resultsOutputInput.value = result.filePaths[0];
        }
      })
      .catch(err => console.error('Error opening folder dialog:', err));
});

const chartTitles= [
    "secretion rate [log(k)]",
    "clearance rate [log(d)]", 
    "pulsatility parameter [f]",
    "pulse duration [log(ton)]",
    "interval between pulses [log(toff)]"
  ];

function createTooltip(element) {
  const tooltipText = element.getAttribute('data-tooltip');
  if (!tooltipText) return;

  const tooltip = document.createElement('div');
  tooltip.className = 'tooltip';
  tooltip.textContent = tooltipText;
  document.body.appendChild(tooltip);

  function positionTooltip(e) {
    const rect = element.getBoundingClientRect();
    const tooltipRect = tooltip.getBoundingClientRect();
    
    let left = rect.left + (rect.width - tooltipRect.width) / 2;
    let top = rect.top - tooltipRect.height - 10;

    // Adjust horizontal position if tooltip goes off-screen
    if (left < 10) left = 10;
    if (left + tooltipRect.width > window.innerWidth - 10) {
      left = window.innerWidth - tooltipRect.width - 10;
    }

    // If tooltip goes off the top, position it below the element
    if (top < 10) {
      top = rect.bottom + 10;
    }

    tooltip.style.left = `${left}px`;
    tooltip.style.top = `${top}px`;
    tooltip.style.opacity = '1';
  }

  element.addEventListener('mouseenter', positionTooltip);
  element.addEventListener('mousemove', positionTooltip);

  element.addEventListener('mouseleave', () => {
    tooltip.style.opacity = '0';
  });
}

// Apply tooltips to all elements with data-tooltip attribute
document.querySelectorAll('[data-tooltip]').forEach(createTooltip);
/** SETTINGS FUNCTIONS */
async function loadSettings() {
  try {
    currentSettings = await window.electron.ipcRenderer.invoke('get-settings');
    if (!currentSettings) {
      throw new Error('No settings found');
    }
    displaySettings();
    // Switch to the Settings tab
    document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
    document.querySelectorAll('.tab-content').forEach(tc => tc.classList.remove('active'));
    document.querySelector('.tab[data-tab="settings"]').classList.add('active');
    document.getElementById('settings').classList.add('active');

  } catch (error) {
    console.error('Error loading settings:', error);
    const container = document.getElementById('settings-container');
    container.innerHTML = `<p>Error loading settings: ${error.message}</p>`;
  }
}

function displaySettings() {
  const container = document.getElementById('settings-container');
  container.innerHTML = '';

  if (!currentSettings) {
    console.error('currentSettings is null or undefined');
    container.innerHTML = '<p>No settings available to display.</p>';
    return;
  }

  for (const [group, settingNames] of Object.entries(settingsToDisplay)) {
    const groupDiv = document.createElement('div');
    groupDiv.className = 'settings-group';
    const groupLabel = groupLabels[group] || group; // Use the descriptive group label if available, otherwise use the group name
    groupDiv.innerHTML = `<h3>${groupLabel}</h3>`;

    settingNames.forEach(name => {
      if (currentSettings[group] && currentSettings[group][name]) {
        const setting = currentSettings[group][name];
        const label = settingLabels[group][name] || name; // Use the descriptive label if available, otherwise use the setting name
        const settingDiv = document.createElement('div');
        settingDiv.className = 'setting-item';
        settingDiv.innerHTML = `
          <label for="${group}___${name}">${label}:</label>
          <input type="text" id="${group}___${name}" value="${setting.value}" data-original-value="${setting.value}">
        `;
        groupDiv.appendChild(settingDiv);
      } else {
        console.warn(`Setting not found: ${group}.${name}`);
      }
    });

    container.appendChild(groupDiv);
  }
  
}

async function showDataFolder() {
  
    const folderPath = document.getElementById('app-data-folder-label').textContent;
    console.log('Attempting to open folder:', folderPath); // Log the path
  
    // Call the openPath function from the exposed electron API
    const a = await await window.electron.ipcRenderer.invoke('open-folder', folderPath);
}

async function showSettingsFolder() {
  
  const folderPath = document.getElementById('app-settings-folder-label').textContent;
  console.log('Attempting to open folder:', folderPath); // Log the path

  // Call the openPath function from the exposed electron API
  const a = await await window.electron.ipcRenderer.invoke('open-folder', folderPath);
}

async function showAppFolder() {
  
  const folderPath = document.getElementById('app-folder-label').textContent;
  console.log('Attempting to open folder:', folderPath); // Log the path

  // Call the openPath function from the exposed electron API
  const a = await await window.electron.ipcRenderer.invoke('open-folder', folderPath);
}


async function saveSettings() {
 
  const inputs = document.querySelectorAll('#settings-container input');
  inputs.forEach(input => {
    const [group, name] = input.id.split('___');
    if (currentSettings[group] && currentSettings[group][name]) {
      currentSettings[group][name].value = input.value;
    } else {
      console.warn(`Setting not found: ${group}.${name}`);
    }
  });

  const success = await window.electron.ipcRenderer.invoke('save-settings', currentSettings);
  if (success) {
    alert('Settings saved successfully!');
    await loadSettings(); // Reload settings after saving
  } else {
    alert('Error saving settings. Please try again.');
  }
}

document.getElementById('save-settings').addEventListener('click', saveSettings);

document.getElementById('settings-container').addEventListener('input', (event) => {
  if (event.target.tagName === 'INPUT') {
    const [group, name] = event.target.id.split('___');
    const originalValue = event.target.getAttribute('data-original-value');
    
  }
});


// Load settings when the page loads
//document.addEventListener('DOMContentLoaded', loadSettings);
//document.querySelector('.tab[data-tab="data"]').classList.add('active');
//document.getElementById('data').classList.add('active');

// Also load settings when the settings tab is clicked
document.querySelector('[data-tab="settings"]').addEventListener('click', loadSettings);

//**SETTINGS FUNCTIONS END**//

document.addEventListener('DOMContentLoaded', () => {
  const loadDataButton = document.getElementById('load-data');
  const loadSettingsButton = document.getElementById('load-settings');
  const startAnalysisButton = document.getElementById('start-analysis');
  const stopAnalysisButton = document.getElementById('stop-analysis');
  const exitAppButton = document.getElementById('exit-app');
  const resultsTab = document.querySelector('.tab[data-tab="results"]');
  const loadResultsButton = document.getElementById('load-results');
  const resultsStringInput = document.getElementById('results-string');
  const resultsNoIdInput = document.getElementById('results-noid');
  const fileList = document.getElementById('file-list');
  window.electron.ipcRenderer.on('log-message', (event, message) => {
    const logPanel = document.getElementById('log-panel');
    const logEntry = document.createElement('div');
    logEntry.textContent = message;
    logEntry.classList.add('log-entry'); // Add class for styling
    logPanel.appendChild(logEntry);
    logPanel.scrollTop = logPanel.scrollHeight; // Auto-scroll to the bottom
  });
  
  document.getElementById('app-data-folder-label').addEventListener('click', showDataFolder);
  document.getElementById('app-settings-folder-label').addEventListener('click', showSettingsFolder);
  document.getElementById('app-folder-label').addEventListener('click', showAppFolder);
  
  document.getElementById('save-settings').addEventListener('click', saveSettings);
  document.getElementById('settings-container').addEventListener('input', (event) => {
    if (event.target.tagName === 'INPUT') {
      const [group, name] = event.target.id.split('___');
      const originalValue = event.target.getAttribute('data-original-value');
      
    }
  });
  resultsTab.addEventListener('click', () => {
    resultsStringInput.value = '';
    fileList.innerHTML = '';
  });
  
  if (loadSettingsButton) {
    loadSettingsButton.addEventListener('click', loadSettings);
  }

   // Load settings when the page loads
  //loadSettings();

  //document.querySelector('.tab[data-tab="data"]').classList.add('active');
  //document.getElementById('data').classList.add('active');

  if (loadResultsButton) {
    loadResultsButton.addEventListener('click', async () => {
      const inputString = resultsStringInput.value.trim() + '_' + resultsNoIdInput.value.trim();
      dataset_id = resultsNoIdInput.value.trim();
      //const res = inputString.split('_');
      //dataset_id = +res[res.length-1];
      const tString = 'lik';
      if (inputString && dataset_id) {
        
        try {
          
          //const files = await window.electron.ipcRenderer.invoke('load-result-files', inputString);
          
          
          // filter files to get the likelihood outputs
          //const likFiles = files.filter((_, index) => index % 3 !== 0 && index % 3 !== 2);
          let res_files = await window.electron.ipcRenderer.invoke('get-all-result-files', inputString, 'lik'); 
          plotLikData(res_files);

          //const trajFiles = files.filter((_, index) => index % 3 !== 0 && index % 3 !== 1);
          
          res_files = await window.electron.ipcRenderer.invoke('get-all-result-files', inputString, 'trajectories'); 
          plotTrajectoryData(res_files, inputString);

          // filter files to get the parameter outputs
          //const paramFiles = files.filter((_, index) => index % 3 !== 1 && index % 3 !== 2);
          res_files = await window.electron.ipcRenderer.invoke('get-all-result-files', inputString, 'parameters'); 
          plotParamData(res_files)
          
          const outputfolder = resultsOutputInput.value.trim();
          const burn_ = document.getElementById('results-burn-in').value;
          const outputPath = await window.electron.ipcRenderer.invoke('write-param-data', res_files, outputfolder, inputString, burn_);
          //writeParamData(paramFiles)

        } catch (error) {
          console.error('Error loading result files:', error);
          alert('An error was encountered while generating the report. Please ensure you are using a valid Job ID and dataset index.');
        }
      } else {
        alert('Please enter a jod identifier and dataset id before loading results.'+inputString + '--'+dataset_id  );
      }
    });
  }

  function displayFileList(files) {
    fileList.innerHTML = '';
    files.forEach(file => {
      const li = document.createElement('li');
      li.textContent = file;
      li.addEventListener('click', () => {
        window.electron.ipcRenderer.send('open-file', file);
      });
      fileList.appendChild(li);
    });
  }

  if (loadDataButton) {
    loadDataButton.addEventListener('click', async () => {
      data = await window.electron.ipcRenderer.invoke('load-data');
      const tableBody = document.getElementById('data-table').querySelector('tbody');
      tableBody.innerHTML = ''; // Clear existing data
      console.log('Loaded data:', data);
      data.forEach((row, index) => {
        const tr = document.createElement('tr');
        const labelTd = document.createElement('td');
        labelTd.textContent = index % 2 === 0 ? `Times ${Math.floor(index / 2) + 1}` : `LH ${Math.floor(index / 2) + 1}`;
        labelTd.classList.add('label-cell');
        tr.appendChild(labelTd);

        Object.values(row).forEach(value => {
          const td = document.createElement('td');
          td.textContent = value;
          tr.appendChild(td);
        });
        tableBody.appendChild(tr);
      });
      
    });
  }
  
  

  document.getElementById('load-settings').addEventListener('click', () => {
    window.electron.ipcRenderer.send('load-settings');
  });


  if (startAnalysisButton) {
    startAnalysisButton.addEventListener('click', () => {
      window.electron.ipcRenderer.send('start-analysis');
      stopAnalysisButton.disabled = false;
    });
  }

  if (stopAnalysisButton) {
    stopAnalysisButton.addEventListener('click',  () => {
      
      window.electron.ipcRenderer.send('stop-analysis');
      stopAnalysisButton.disabled = true;
      
    });
  }

  

  if (exitAppButton) {
    exitAppButton.addEventListener('click', () => {
      window.electron.ipcRenderer.send('exit-app');
    });
  }
  


  // Tab switching logic
  document.querySelectorAll('.tab').forEach(tab => {
    tab.addEventListener('click', () => {
      document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
      document.querySelectorAll('.tab-content').forEach(tc => tc.classList.remove('active'));

      tab.classList.add('active');
      document.getElementById(tab.getAttribute('data-tab')).classList.add('active');
    });
  });

 
});

async function plotLikData(files) {
  try {
    const burn_ = document.getElementById('results-burn-in').value;
    const chartData = await window.electron.ipcRenderer.invoke('pool-and-plot-data', files, burn_);
    const ctx = document.getElementById('pooled-chart').getContext('2d');
    
    console.log(chartData.labels.map((y, index) => ({x: parseInt(index), y: parseFloat(y)})));
    // Destroy the existing chart if it exists
    if (currentChart) {
      currentChart.destroy();
    }
    
    // Create a new chart and store the reference
    currentChart = new Chart(ctx, {
      type: 'line',
      data: {
      datasets: [
        {
          data: chartData.labels.map((y, index) => ({x: parseInt(index), y: parseFloat(y)})),
          borderWidth: 1,
          borderColor: 'rgba(0, 0, 0, 1)',
          backgroundColor: 'rgba(0, 0, 0, 0.2)',
          pointRadius: 0,
          fill: false,
          showLine: true
        }]},
        options: {
          
            responsive: false,
            maintainAspectRatio: false,
            plugins: {
              legend: { display: false }
            },
            scales: {
              x: {
                type: 'linear',
                display: true,
                title: {
                  display: true,
                  text: 'mcmc iter',
                  color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
                },
                ticks: { 
                  display: true,
                  color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
                 }
                
              },
              y: {
                type: 'linear',
                display: true,
                title: {
                  display: true,
                  text: 'log-likelihood',
                  color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
                },
                grid: {
                  color: 'rgba(0, 0, 0, 0.1)' // Light grid lines// 'rgba(255, 255, 255, 0.1)' // Light grid lines
                },
                ticks: {
                  color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
                }
                
              }
            }
            
          }
      
    });
  } catch (error) {
    console.error('Error plotting likelihood data:', error);
    
  }
}
/* plot trajectory data */
async function plotTrajectoryData(files, inputString) {
  
  try {
    const burn_ = document.getElementById('results-burn-in').value;
    const trajectoryData = await window.electron.ipcRenderer.invoke('load-and-process-trajectory-data', files,burn_);
    

    
    // Process the data
    const timeData = trajectoryData[0].timeData; // Assuming all files have the same time points
    
    // Concatenate data from all files
    const allGenData = trajectoryData.flatMap(file => file.genData);
    const allLHData = trajectoryData.flatMap(file => file.LHData);

    
    // Calculate averages
    const avgGenData = calculateAverage(allGenData);
    const avgLHData = calculateAverage(allLHData);

    
    
    // Find indices where avgGenData crosses 0.5 from below
    const threshold = getThresholdValue();
    
    if (threshold !== null) {
      console.log('---------------->> Valid threshold:', threshold);  
      
      const crossingIndices = findCrossingIndicesFromBelow(avgGenData, threshold);
      // Find corresponding time values
      const crossingTimes = crossingIndices.map(index => timeData[0][index]);
      let LHvals = []
      for (let i = 0; i < crossingIndices.length; i++) {
        
        
        let max = 0;
        let min = 100;
        let startIndex = crossingIndices[i] - 10;
        let endIndex = crossingIndices[i] + 15;
        if (endIndex >= timeData[0].length) {
          endIndex = timeData[0].length-1;  
        }
        // Loop through the specified range
        for (let iii = startIndex; iii <= endIndex; iii++) {
          
          if (avgLHData[iii] > max) {
            max = avgLHData[iii] ; // Update max if a larger value is found
          }
          if (avgLHData[iii] < min) {
            min = avgLHData[iii] ; // Update max if a larger value is found
          }
        }
        // Update max if a larger value is found
        
        LHvals.push(max -min);
      }


      let i = 0;
      let timesData2 = [];
      let LHData = [];
     ;

      if (dataset_id !== null && data !== null && data.length > 0) {
        const requestedRaw = String(dataset_id ?? '').trim();
        // filter rows whose first column matches dataset_id (robust compare)
        const matchedRows = data.filter(row => {
          const firstRaw = row[0];
          return firstRaw === requestedRaw;
        });

        if (matchedRows.length >= 2) {
          timesData2 = matchedRows[0];
          LHData    = matchedRows[1];
        }
        
        
 
      }

      
      // Plot the data
      plotAverageData('gen-chart', timeData[0], avgGenData, 'Hypothalamic drive',[],[], []);
      plotAverageData('lh-chart', timeData[0], avgLHData, 'LH', crossingTimes, timesData2, LHData);

      // Create a new 2D array
      const combinedArray = [];

      // Loop through the arrays and combine them
      for (let i = 0; i < crossingTimes.length; i++) {
          combinedArray.push([crossingTimes[i], LHvals[i]]); // Create a new inner array for each pair
      }  
      // Write crossing times to file
      await writeCrossingTimesToFile(combinedArray, inputString);
    }

  } catch (error) {
    console.error('Error plotting trajectory ', error);
    
  }
}

// Function to get and validate the threshold value
function getThresholdValue() {
  // Get the input value from the results-threshold input field
  const thresholdInput = document.getElementById('results-threshold').value;

  // Convert the input value to a number
  const thresholdValue = parseFloat(thresholdInput);

  // Validate the number
  if (isNaN(thresholdValue)) {
    alert('Please enter a valid threshold.'); // Alert if not a number
    return null; // Return null or handle the error as needed
  }

  if (thresholdValue < 0 || thresholdValue > 1) {
    alert('Please enter a threshold between 0 and 1.'); // Alert if out of range
    return null; // Return null or handle the error as needed
  }

  return thresholdValue; // Return the valid threshold value
}

async function writeCrossingTimesToFile(crossingTimes, inputString) {
  const fileName = `${inputString}_pulse_times.csv`;
  const csvContent = crossingTimes.join('\n');
  const outputfolder = resultsOutputInput.value.trim();
  try {
    await window.electron.ipcRenderer.invoke('write-file', fileName, outputfolder, csvContent);
    
  } catch (error) {
    console.error('Error writing crossing times to file:', error);
  }
}
function findCrossingIndicesFromBelow(data, threshold) {
  const crossingIndices = [];
  for (let i = 0; i < data.length; i++) {
    if (i==0) {
      if (data[i] >= threshold) {
        crossingIndices.push(i);
      }
    }
    else if (data[i-1] < threshold && data[i] >= threshold) {
      crossingIndices.push(i);
    }
  }
  return crossingIndices;
}

function calculateAverage(data) {
  // Initialize an array to store the sum of each column
  const sum = new Array(data[0].length).fill(0);

  // Sum up all values for each column
  data.forEach(row => {
    row.forEach((value, index) => {
      sum[index] += parseFloat(value);
    });
  });

  // Divide each sum by the number of rows to get the average
  return sum.map(columnSum => columnSum / data.length);
}

function plotAverageData(canvasId, timeData, modelData, label, crossingTimes, timesData2, LHData) {
  const ctx = document.getElementById(canvasId).getContext('2d');
  
  // Destroy existing chart if it exists
  if (window.trajectoryCharts && window.trajectoryCharts[canvasId]) {
    window.trajectoryCharts[canvasId].destroy();
  }
  
  
 //const timesData2 = data[(i-1)*2];
  //const LHData = data[(i-1)*2+1];
  const maxY = Math.max(...modelData, ...Object.values(LHData));
  
  // Create new chart
  window.trajectoryCharts = window.trajectoryCharts || {};
  window.trajectoryCharts[canvasId] = new Chart(ctx, {
    type: 'line',
    data: {
      //labels: timesData,
      datasets: [
      {
        label: 'Model',
        data: modelData.map((y, index) => ({x: parseFloat(timeData[index]), y: parseFloat(y)})),
        borderColor: 'rgba(0, 0, 0, 1)',
        backgroundColor: 'rgba(0, 0, 0, 0.2)',
        borderWidth: 2,
        pointRadius: 0,
        fill: false,
        order: 0
      },
      {
        label: 'Data',
        data: Object.values(LHData).map((y, index) => ({x: parseFloat(timesData2[index]), y: parseFloat(y)})),
        borderColor: 'rgba(0, 0, 0, 1)',
        backgroundColor: 'rgba(0, 0, 0, 0.2)',
        borderWidth: 1,
        pointRadius: 2,
        fill: false,
        showLine: true,
        order: 1
      },
      {
        label: 'Pulses',
        data: crossingTimes.map((time,index) => ({x: parseFloat(time), y: parseFloat(maxY)})),
        borderColor: 'rgba(255, 0, 0, 1)',
        backgroundColor: 'rgba(255, 0, 0, 1)',
        borderWidth: 1,
        pointRadius: 3,
        pointStyle: 'star',
        showLine: false,
        order: 2
      }
    ]
    },
    options: {
      responsive: false,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          display: false // Hide the legend if you don't need it
        }
      },  
      scales: {
        x: {
          type: 'linear',
          position: 'bottom',
          title: {
            display: true,
            text: 'Time (min)',
            color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
          },
          grid: {
            color: 'rgba(0, 0, 0, 0.1)' // Light grid lines
          },
          ticks: {
            color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
          }
        },
        y: {
          type: 'linear',
          position: 'left',
          title: {
            display: true,
            text: label,
            color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
          },
          grid: {
            color: 'rgba(0, 0, 0, 0.1)' // Light grid lines
          },
          ticks: {
            color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
          }
        }
      }
    }
  });

}
/* end plot trajectory data  */

// Download the plot as an image
document.querySelectorAll('.download-plot').forEach(button => {
  button.addEventListener('click', (event) => {
    const canvasId = event.target.getAttribute('data-canvas-id');
    const canvas = document.getElementById(canvasId);
    const link = document.createElement('a');
    
    // Get values from input fields
    const jobId = document.getElementById('results-string').value; // Get value from results-string input
    const datasetIndex = document.getElementById('results-noid').value; // Get value from results-noid input

    // Create the filename using the input values
    const filename = `${canvasId}_${jobId}_${datasetIndex}.png`; // Example: plot_123_0.png


    link.href = canvas.toDataURL('image/png'); // Convert canvas to data URL
    link.download = filename; // Set the default filename
    link.click(); // Trigger the download
  });
});

async function plotParamData(files) {

  try {
    const burn_ = document.getElementById('results-burn-in').value;
    const chartsData = await window.electron.ipcRenderer.invoke('pool-and-plot-parameter-data', files, burn_);
    //console.log('Received chartsData:', chartsData);

    // Clear existing charts
    currentCharts.forEach(chart => chart.destroy());
    currentCharts = [];

    const chartContainer = document.getElementById('charts-container');
    chartContainer.innerHTML = '';

   

    // Create new charts
    chartsData.forEach((chartData, index) => {
      
      const wrapper = document.createElement('div');
      wrapper.className = 'plot-container';
      
      const canvas = document.createElement('canvas');
      canvas.id = `chart-${index}`;
      

      const button = document.createElement('button');
      button.className = 'download-plot';
      button.setAttribute('data-canvas-id', canvas.id);
      button.style.display = 'block'; // Initially hidden
      button.innerHTML = '<i class="fas fa-download"></i>'; // Font Awesome icon

      

      wrapper.appendChild(canvas);
      wrapper.appendChild(button)
      chartContainer.appendChild(wrapper);

      const ctx = canvas.getContext('2d');
      const newChart = new Chart(ctx, {
        type: 'line',
        data: {
          labels: chartData.labels,
          datasets: [{
            ...chartData.datasets[0],
            borderColor: 'rgba(0, 0, 0, 1)',
            backgroundColor: 'rgba(0, 0, 0, 0.8)',
            borderWidth: 1,
            pointRadius: 0,
            fill: false
          }]
        },
        options: {
          
          responsive: false,
          maintainAspectRatio: true,
          plugins: {
            title: {
              display: true,
              text: chartTitles[index] || `Parameter ${index + 1}`,
              font: { size: 12 },
              color: 'rgba(0, 0, 0, 0.8)' 
            },
            legend: { display: false }
          },
          scales: {
            x: {
              type: 'linear',
              display: true,
              ticks: { 
                display: true, 
                color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
               },
              title: {
                display: true,
                text: 'mcmc iter',
                color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
              }
            },
            y: {
              display: true,
              title: { display: false },
              grid: {
                color: 'rgba(0, 0, 0, 0.1)' // Light grid lines
              },
              ticks: {
                color: 'rgba(0, 0, 0, 0.8)' // Light text color for ticks
              }
              
            }
          }
          
        }
      });
      currentCharts.push(newChart);
      // Add event listener for the download button
      button.addEventListener('click', () => {
        const link = document.createElement('a');
        // Get values from input fields
        const jobId = document.getElementById('results-string').value; // Get value from results-string input
        const datasetIndex = document.getElementById('results-noid').value; // Get value from results-noid input

        // Create the filename using the input values
        const filename = `chart-${index}_${jobId}_${datasetIndex}.png`; // Example: plot_123_0.png

        link.href = canvas.toDataURL('image/png'); // Convert canvas to data URL
        link.download = filename; // Set the default filename
        link.click(); // Trigger the download
      });

    });

    console.log(`Created ${currentCharts.length} charts`);
  } catch (error) {
    console.error('Error plotting parameter data:', error);
    
  }
  
}


const version = window.electron.getAppVersion();
window.electron.getAppVersion().then(version => {
  document.title = `HormoneBayes Analysis Tool (v${version})`;
});




// Replace the click handler with mousedown and a guard to avoid rebuilding repeatedly
document.getElementById('results-string').addEventListener('pointerdown', async (e) => {
  const dropdown = document.getElementById('results-string');
  const folder = document.getElementById('app-folder-label');

  try {
    const files = await window.electron.ipcRenderer.invoke('read-outputs', folder.textContent);
    const fileKey = files.join('\n');

    // If file list unchanged and dropdown already populated, just open it
    if (dropdown.dataset._files === fileKey && dropdown.options.length > 1) {
      dropdown.focus();
      return;
    }

    // cache file list
    dropdown.dataset._files = fileKey;

    // build options before the native menu opens
    dropdown.innerHTML = '';
    const placeholderOption = document.createElement('option');
    placeholderOption.value = '';
    placeholderOption.textContent = 'Select Job ID';
    placeholderOption.disabled = true;
    placeholderOption.selected = true;
    dropdown.appendChild(placeholderOption);

    const uniqueNames = new Map();
    files.forEach(file => {
      if (!file.endsWith('.csv')) return;
      const parts = file.split('_');
      const displayName = parts.length > 4 ? parts.slice(0, parts.length - 4).join('_') : file;
      const datasetIndex = parts[Math.max(0, parts.length - 4)];
      const list = uniqueNames.get(displayName) || [];
      list.push(datasetIndex);
      uniqueNames.set(displayName, list);
    });

    uniqueNames.forEach((value, key) => {
      const se = new Set(value);
      const option = document.createElement('option');
      option.value = key;
      option.textContent = `${key} (${se.size} datasets)`;
      dropdown.appendChild(option);
    });

    // ensure the select receives focus so native menu opens
    dropdown.focus();
  } catch (err) {
    console.error('Error reading outputs folder:', err);
  }
});


document.getElementById('results-string').addEventListener('change', (event) => {
  const selectedValue = event.target.value; // Get the selected value
  

  // You can perform additional actions based on the selection
  if (selectedValue) {
      // Do something with the selected value
      const val = new Set(uniqueNames.get(selectedValue));
      const d_input = document.getElementById('results-noid');
      let text = "";
      val.forEach (function(value) {
        text += `${value}, `;
      })
      text = text.slice(0,-2);
      
      d_input.placeholder = `Valid indices: ${text}.`;
      //alert(`You selected: ${selectedValue}`);
  }
});



