const { app, BrowserWindow, ipcMain, dialog,shell } = require('electron');
const path = require('path');
const fs = require('fs');
const csv = require('csv-parser');
const { exec, spawn } = require('child_process');
const os = require('os');
const fsp = require('fs').promises;

let mainWindow;
let dataWindow;
let csvData = [];
const hb_container = `hormonebayes:${app.getVersion()}`;

const burn_prc = 10; // Change this to the percentage you want to remove

const parameter_number = 5;

const options = {
  env: {
    PATH: `${process.env.PATH}:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin`
  }
};

// Add this near the top of your main.js file
const chartTitles = [
  "secretion rate (log k)",
  "clearance rate (log k)", 
  "pulsatility parameter (f)",
  "pulse duration (log k)",
  "interval between pulses (log k)"
];

function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1000,
    height: 800,
    title: "HormoneBayes Analysis Tool",
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      enableRemoteModule: false,
      nodeIntegration: false
    }
  });

  mainWindow.loadFile('index.html');
  
  // Send the version to the renderer process
  mainWindow.webContents.on('did-finish-load', () => {
    mainWindow.webContents.send('app-version', app.getVersion());
  });

  // Send the app path to the renderer process
  mainWindow.webContents.on('did-finish-load', () => {
    const appPath = path.join(app.getAppPath(), 'out');
    mainWindow.webContents.send('app-out-path', appPath);
  });

  // Send the app path to the renderer process
  mainWindow.webContents.on('did-finish-load', () => {
    const appPath = path.join(app.getAppPath(), 'data');
    mainWindow.webContents.send('app-data-path', appPath);
  });

  // Open the DevTools.
  //mainWindow.webContents.openDevTools();

  // Handle window close event
  mainWindow.on('closed', () => {
    mainWindow = null;
  });
}

ipcMain.handle('get-app-version', () => app.getVersion());

function createDataWindow() {
  dataWindow = new BrowserWindow({
    width: 1000,
    height: 800,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      enableRemoteModule: false,
      nodeIntegration: false
    }
  });

  dataWindow.loadFile('data.html');

  // Handle window close event
  dataWindow.on('closed', () => {
    dataWindow = null;
  });

  // Send data to the data window
  dataWindow.webContents.on('did-finish-load', () => {
    dataWindow.webContents.send('send-data', csvData);
    dataWindow.webContents.send('app-version', app.getVersion());
  });
}

function logMessage(message) {
  const timestamp = new Date().toLocaleString();
  const formattedMessage = `> ${message} [${timestamp}]`;
  console.log(formattedMessage);
  
  mainWindow.webContents.send('log-message', formattedMessage);
  
}

function checkDockerDaemon() {
  exec('docker info',options, (error, stdout, stderr) => {
    if (error) {
      logMessage('Docker daemon is not running');
      if (mainWindow) {
        dialog.showMessageBox(mainWindow, {
          type: 'warning',
          buttons: ['Start Docker', 'Cancel'],
          defaultId: 0,
          title: 'Docker Daemon',
          message: 'Docker daemon is not running. Would you like to start it?',
        }).then(result => {
          if (result.response === 0) {
            startDockerDaemon();
          }
        });
      }
    } else {
      logMessage('Docker daemon is running');
    }
  });
}

function startDockerDaemon() {
  const platform = os.platform();
  let command;

  if (platform === 'win32') {
    command = 'powershell.exe -Command "Start-Process \'C:\\Program Files\\Docker\\Docker\\Docker Desktop.exe\'"';
  } else if (platform === 'darwin') {
    command = 'open --background -a Docker';
  } else if (platform === 'linux') {
    command = 'sudo service docker start';
  } else {
    logMessage(`Unsupported platform: ${platform}`);
    return;
  }

  

  exec(command, options, (error, stdout, stderr) => {
    if (error) {
      logMessage('Failed to start Docker. Try starting manually.');
    } else {
      logMessage('Docker daemon started successfully');
    }
  });
}

app.whenReady().then(() => {
  createWindow();
  checkDockerDaemon();

  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) {
      createWindow();
    }
  });
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

ipcMain.handle('load-data', async () => {
  


  const result = await dialog.showOpenDialog(mainWindow, {
    properties: ['openFile'],
    filters: [{ name: 'CSV Files', extensions: ['csv'] }]
  });

  if (result.canceled) {
    return [];
  }

  const filePath = result.filePaths[0];
  const appPath = app.getAppPath();
  const dataDir = path.join(appPath, 'data');
  const destPath = path.join(dataDir, 'data.csv');
  
  // Ensure the data directory exists
  if (!fs.existsSync(dataDir)) {
    fs.mkdirSync(dataDir);
  }


  // Copy the CSV file to the data directory
  fs.copyFileSync(filePath, destPath);
  

  csvData = [];

  return new Promise((resolve, reject) => {
    fs.createReadStream(filePath)
      .pipe(csv({ headers: false })) // Treat all rows as data
      .on('data', (row) => {
        csvData.push(row);
      })
      .on('end', () => {
        logMessage(`Data loaded successfully and saved in ${destPath}`);
        if (mainWindow && !mainWindow.isDestroyed()) {
          mainWindow.webContents.send('send-data', csvData); // Send data to the renderer process
        }
        resolve(csvData);
      })
      .on('error', (error) => {
        reject(error);
      });
  });
});

ipcMain.on('see-data', () => {
  if (!dataWindow) {
    createDataWindow();
  } else {
    dataWindow.focus();
    // Send data to the data window if it already exists
    if (dataWindow && !dataWindow.isDestroyed()) {
      dataWindow.webContents.send('send-data', csvData);
    }
  }
});

ipcMain.handle('get-data', () => {
  return csvData;
});

ipcMain.on('close-window', () => {
  if (mainWindow) {
    mainWindow.close();
  }
});

ipcMain.on('log-message', (event, message) => {
  logMessage(message);
});

// Handle new IPC messages
ipcMain.on('load-settings', async () => {
  
  const result = await dialog.showOpenDialog(mainWindow, {
    properties: ['openFile'],
    filters: [{ name: 'Text Files', extensions: ['txt'] }]
  });

  if (result.canceled) {
    return;
  }

  const filePath = result.filePaths[0];
  const appPath = app.getAppPath();
  const dataDir = path.join(appPath, 'data');
  const destPath = path.join(dataDir, 'params.txt');

  // Ensure the data directory exists
  if (!fs.existsSync(dataDir)) {
    fs.mkdirSync(dataDir);
  }

  // Copy the TXT file to the data directory and rename it to params.txt
  fs.copyFileSync(filePath, destPath);
  logMessage(`Settings loaded successfully and saved in ${destPath}`);
});

ipcMain.on('view-settings', () => {
  // Add your logic for viewing settings here
});

ipcMain.on('start-analysis', () => {
  logMessage('Building docker container ... ');
  const outPath = path.join(app.getAppPath(), 'out');
  
  const buildProcess = spawn('docker', ['build',  '-t', hb_container, app.getAppPath()],options);

  buildProcess.stdout.on('data', (data) => {
    logMessage(`${data}`);
  });

  buildProcess.stderr.on('data', (data) => {
    logMessage(`Build stderr: ${data}`);
  });

  buildProcess.on('close', (code) => {
    if (code !== 0) {
      logMessage(`Build process exited with code ${code}`);
      return;
    }
    logMessage('Docker image built successfully');

    const runProcess = spawn('docker', ['run', '--rm', '-v', `${outPath}:/usr/src/myapp/out`, '--name', 'my_analysis_container', hb_container],options);

    runProcess.stdout.on('data', (data) => {
      logMessage(`${data}`);
    });

    runProcess.stderr.on('data', (data) => {
      logMessage(`Run stderr: ${data}`);
    });

    runProcess.on('close', (code) => {
      if (code !== 0) {
        logMessage(`Run process exited with code ${code}`);
        return;
      }
      logMessage('Analysis finished successfully');
    });
  });
});

ipcMain.on('stop-analysis', async () => {
  logMessage('Attempting to stop analysis ... ');
  const killCommand = 'docker kill my_analysis_container';
  
  
  const code1 = exec(killCommand, options,  (killError, killStdout, killStderr) => {
    if (killError) {
      logMessage(`Failed to stop Docker container: ${killStderr}`);
      return -1;
    }
    logMessage('Analysis stopped successfully');
    return 0;
    
  });
  console.log(`code is ${code1}`);
  return code1;
});

ipcMain.on('exit-app',  () => {
  logMessage('Checking if anything is still running ... ');
  const killCommand = 'docker kill my_analysis_container';
  

  exec(killCommand, options, (killError, killStdout, killStderr) => {
    if (killError) {
      logMessage(`Nothing to stop`);
      return;
    }
    logMessage('Analysis stopped successfully');
  });
  app.quit();
});


ipcMain.handle('load-result-files', async (event, inputString) => {
  const baseDir = path.join(__dirname, 'out');
  const files = [];
  let c1 = 0;
  while (true) {
    const filePatterns = 
     [`${inputString}_chain_${c1}_parameters.csv`,
      `${inputString}_chain_${c1}_lik.csv`,
      `${inputString}_chain_${c1}_trajectories.csv`
    ];

    const existingFiles = await Promise.all(
      filePatterns.map(async (pattern) => {
        const filePath = path.join(baseDir, pattern);
        
        try {
          await fsp.access(filePath);
          return filePath;
        } catch {
          return null;
        }
      })
    );
    //onsole.log( existingFiles)
    const validFiles = existingFiles.filter(Boolean);
    if (validFiles.length === 0) break;

    files.push(...validFiles);
    c1++;
  }
  
  return files;
});


ipcMain.handle('get-all-result-files', async (event, inputString, ending) => {
  const baseDir = path.join(__dirname, 'out');
  const allFiles = await fs.promises.readdir(baseDir); // Read the directory
  const regex = new RegExp(`${inputString}_chain_[0-9]+_${ending}.csv`, 'g'); // Create a regex from the pattern
  const files = [];
  // Count occurrences of the pattern in the filenames
  
  allFiles.forEach(file => {
      if (file.match(regex)) {
        files.push(path.join(baseDir, file));
      }
    
  });
  
  return files;
});







ipcMain.on('open-file', (event, filePath) => {
  const command = process.platform === 'win32' ? 'start' : process.platform === 'darwin' ? 'open' : 'xdg-open';
  exec(`${command} "${filePath}"`, (error) => {
    if (error) {
      console.error(`Error opening file: ${error}`);
    }
  });
});



async function poolAndPlotData(filenames, burn_) {
  try {
    const labels = [];
    const datasets = [];


    for (const filename of filenames) {
      // Read the file content
      const fileContent = await fsp.readFile(filename, 'utf8');
      
      // Split the content into lines and parse each line as a float
      const allValues = fileContent.trim().split('\n').map(line => parseFloat(line.trim()));
      
      // Calculate how many values to discard
      const discardCount = Math.floor(allValues.length * (burn_ / 100));
      
      // Slice the array to remove the first a% of values
      const keptValues = allValues.slice(discardCount);
      
      // Add the kept values to the labels array
      labels.push(...keptValues);
    }

    return { labels, datasets};
  } catch (error) {
    console.error('Error processing files:', error);
    throw error;
  }
  
}

// Add this to your existing IPC handlers
ipcMain.handle('pool-and-plot-data', async (event, filenames, burn_) => {
  try {
    
    const chartData = await poolAndPlotData(filenames, burn_);
    return chartData;
  } catch (error) {
    console.error('Error pooling and plotting data:', error);
    throw error;
  }
});


async function readCSVFile(filePath) {
  return new Promise((resolve, reject) => {
    const results = [];
    fs.createReadStream(filePath)
      .pipe(csv({ headers: false }))
      .on('data', (data) => results.push(data))
      .on('end', () => resolve(results))
      .on('error', reject);
  });
}

async function poolAndPlotParameterData(filenames, burn_) {
  const allData = [];
  
  for (const filename of filenames) {
    try {
      console.log(`Reading parameter file: ${filename}`);
      const fileData = await readCSVFile(filename);
      allData.push(fileData);
    } catch (error) {
      console.error(`Error reading file ${filename}:`, error);
      throw error;
    }
  }

  console.log(`Reading parameter ${allData.length}`)
  // Process data into chart format
  const chartsData = processDataIntoCharts(allData, burn_);
  
  //console.log('Processed charts data:', chartsData);

  return chartsData;
}

function processDataIntoCharts(allData, burn_) {
  //console.log('Processing data into charts. Input data:', allData);
  
  if (!Array.isArray(allData) || allData.length === 0) {
    console.error('Invalid input data');
    return [];
  }

  // Initialize chartsData with 5 datasets (one for each parameter)
  const chartsData = new Array(parameter_number).fill(null).map(() => ({
    labels: [],
    datasets: [{ data: [], label: 'Parameter' }],
    means:[],
    vars:[],
    p_mean:0
  }));


  allData.forEach((fileData, fileIndex) => {
    //console.log(`Processing file ${fileIndex}:`, fileData);
    if (!Array.isArray(fileData) || fileData.length !== parameter_number) {
      console.warn(`Skipping invalid file data at index ${fileIndex}`);
      return;
    }

    fileData.forEach((row, rowIndex) => {
      //console.log(`Processing row ${rowIndex}:`, row);
      if (typeof row !== 'object' || row === null) {
        console.warn(`Skipping invalid row at index ${rowIndex} in file ${fileIndex}`);
        return;
      }

      
      const aa = chartsData[rowIndex].labels.length
      const bb = Object.entries(row).length
      const a = burn_; // Change this to the percentage you want to remove

      Object.entries(row).forEach(([key, value]) => {
        const parsedValue = parseFloat(value);
        if (!isNaN(parsedValue)) {
          
          if(key >= bb*(a/100)) {
            chartsData[rowIndex].labels.push(parseInt(key) + aa - bb*(a/100) );
          
            chartsData[rowIndex].datasets[0].data.push(parsedValue);
          }
          
        } else {
          console.warn(`Skipping invalid value: ${value} at row ${rowIndex}, key ${key} in file ${fileIndex}`);
        }
      });

      const ad = bb*(100-a)/100
      chartsData[rowIndex].vars.push( calculateVariance(chartsData[rowIndex].datasets[0].data.slice(-ad)) )
      chartsData[rowIndex].means.push( calculateMean(chartsData[rowIndex].datasets[0].data.slice(-ad)) )

    });
  });

  chartsData.forEach((chartData, index) => {

    chartsData[index].p_mean = calculateMean(chartData.means);
    const N = chartData.datasets[0].data.length / chartData.means.length
    const M = chartData.means.length;
    
    const squaredDifferences = chartsData[index].means.map(val => Math.pow(val - chartsData[index].p_mean, 2)); // Step 2: Calculate squared differences
    const B = squaredDifferences.reduce((acc, val) => acc + val, 0) * N / (M-1); // Step 3: Calculate variance
    const W = calculateMean(chartData.vars)
    const V = (N-1)/N*W + (M+1)/M/N*B
    chartsData[index].R_hat = Math.sqrt(V/W)

    console.log(`Rhat for parameter ${index} is: ${chartsData[index].R_hat}, N=${N}, M=${M}` )
    console.log(`means ${chartData.means} vars: ${chartData.vars}` )
  });
  
   
  //console.log('Processed charts data:', chartsData);
  return chartsData;
}

//intra_chain_means = squeeze(nanmean(dd));
//intra_chain_variances = squeeze(nanvar(dd));
//pooled_mean = mean(intra_chain_means);

//B = N/(M-1) * sum((intra_chain_means-pooled_mean).^2);
//W = mean(intra_chain_variances);
//V = (N-1)/N*W + (M+1)/M/N*B;

//R_hat = sqrt(V/W);

function calculateVariance(data) {
  const mean = data.reduce((acc, val) => acc + val, 0) / data.length; // Step 1: Calculate mean
  const squaredDifferences = data.map(val => Math.pow(val - mean, 2)); // Step 2: Calculate squared differences
  const variance = squaredDifferences.reduce((acc, val) => acc + val, 0) / data.length; // Step 3: Calculate variance
  return variance;
}



ipcMain.handle('pool-and-plot-parameter-data', async (event, filenames, burn_) => {
  try {
    const chartsData = await poolAndPlotParameterData(filenames, burn_);
    return chartsData;
  } catch (error) {
    console.error('Error pooling and plotting parameter data:', error);
    throw error;
  }
});


function calculateHPD(data, alpha) {
  // Sort the data in ascending order
  const sortedData = data.slice().sort((a, b) => a - b);
  const n = sortedData.length;
  
  // Calculate the number of data points to include in the HPD interval
  const intervalSize = Math.floor(n * (1 - alpha));
  
  let minRange = Infinity;
  let hpdMin, hpdMax;

  // Slide a window of size intervalSize through the sorted data
  for (let i = 0; i <= n - intervalSize; i++) {
    const currentMin = sortedData[i];
    const currentMax = sortedData[i + intervalSize - 1];
    const currentRange = currentMax - currentMin;

    // If this range is smaller than the current minimum range, update the HPD interval
    if (currentRange < minRange) {
      minRange = currentRange;
      hpdMin = currentMin;
      hpdMax = currentMax;
    }
  }

  return [hpdMin, hpdMax];
}


// Add this to your existing IPC handlers
ipcMain.handle('calculate-hpd', (event, data, alpha) => {
  return calculateHPD(data, alpha);
});


async function writeParamData(files, outputfolder, outputPrefix, burn_) {
  console.log('Starting writeParamData');
  try {
    const chartsData = await poolAndPlotParameterData(files, burn_);
    console.log('Received chartsData:', chartsData);

    if (!Array.isArray(chartsData) || chartsData.length === 0 || !chartsData[0].datasets) {
      throw new Error('No valid chart data received');
    }

    let outputContent = '';

    const parameterNames = [
      "secretion rate [log(k)]",
      "clearance rate [log(d)]", 
      "pulsatility parameter [f]",
      "pulse duration [log(ton)]",
      "interval between pulses [log(toff)]"
    ];
    // Remove initial a% of the data
    

    chartsData.forEach((chartData, index) => {
      if (!chartData.datasets || !chartData.datasets[0] || !Array.isArray(chartData.datasets[0].data)) {
        console.warn(`Skipping invalid chart data for parameter ${index + 1}`);
        return;
      }

      let parameterData = chartData.datasets[0].data;
      const removeCount = Math.floor(parameterData.length * (burn_ / 100));
      //parameterData = parameterData.slice(removeCount);

      // Apply transformation to the third dataset (index 2)
      if (index === 2) {
        parameterData = parameterData.map(x => 1 / (1 + Math.exp(-x)));
      }

      const parameterName = parameterNames[index] || `Parameter ${index + 1}`;

      // Calculate statistics
      const median = calculateMedian(parameterData);
      const mean = calculateMean(parameterData);
      const [hpdMin, hpdMax] = calculateHPD(parameterData, 0.05); // 95% HPD

      // Append to output content
      outputContent += `${parameterName}:\n`;
      outputContent += `  Median: ${median.toFixed(4)}\n`;
      outputContent += `  Mean: ${mean.toFixed(4)}\n`;
      outputContent += `  95% HPD: [${hpdMin.toFixed(4)}, ${hpdMax.toFixed(4)}]\n`;
      outputContent += `  Rhat: ${chartData.R_hat}\n\n`;
    });

    if (outputContent === '') {
      throw new Error('No valid data to write');
    }

    
    
    
    const t1 = chartsData[3].datasets[0].data.map(x => 10**(x));
    const t2 = chartsData[4].datasets[0].data.map(x => 10**(x));
    const f_param_values = chartsData[2].datasets[0].data.map(x => 1.0/(1.0 + Math.exp(-x)) );
    const dd_values = chartsData[1].datasets[0].data.map(x => 10**(x));
    const a1_values = chartsData[0].datasets[0].data.map(x => 10**(x));

    const dd_mean = calculateMean(dd_values);
    const dd_median = calculateMedian(dd_values);
    const [dd_hpd_min, dd_hpd_max] = calculateHPD(dd_values, 0.05);  

    // Create a new array to hold the results
    const freq_res = [];

    // Loop through the arrays and add elements
    for (let i = 0; i < t1.length; i++) {
      freq_res.push(60/(t1[i] + t2[i]));
    }

    const freq_median = calculateMedian(freq_res);
    const freq_mean = calculateMean(freq_res);
    const [freq_hpd_min, freq_hpd_max] = calculateHPD(freq_res, 0.05); // 95% HPD

    
    //const amplitude = (dd)*(1.0*f_param1 + 1/2.0*(1.0-f_param1)  )/20.0 + 
    //             Math.exp(-a1*tau_on)*(dd)/(20.0/a1-1.0)*(1/20.0-(1.0*f_param1 + 1/2.0*(1.0-f_param1)  )/a1) + 
    //Math.exp(-20.0*tau_on)*(dd)/(20.0/a1-1.0)*((1.0*f_param1 + 1/2.0*(1.0-f_param1)  )/20.0-1/20.0)
    const amplitude = [];
    let   a_val;
    // Loop through the arrays and add elements
    for (let i = 0; i < t1.length; i++) {
      a_val  = (a1_values[i])*(1.0*f_param_values[i] + 1/2.0*(1.0-f_param_values[i])  )/20.0 + 
                    Math.exp(-dd_values[i]*t1[i])*(a1_values[i])/(20.0/dd_values[i]-1.0)*(1/20.0-(1.0*f_param_values[i] + 1/2.0*(1.0-f_param_values[1])  )/dd_values[i]) + 
                    Math.exp(-20.0*t1[i])*(a1_values[i])/(20.0/dd_values[i]-1.0)*((1.0*f_param_values[1] + 1/2.0*(1.0-f_param_values[1])  )/20.0-1/20.0);
      a_val = (a1_values[i]/20./dd_values[i]  ) *(1 -Math.exp(-t1[i]*dd_values[i])) *f_param_values[i] ;
      amplitude.push(a_val);
    }

    const amplitude_median = calculateMedian(amplitude);
    const amplitude_mean = calculateMean(amplitude);
    const [amplitude_hpd_min, amplitude_hpd_max] = calculateHPD(amplitude, 0.05); // 95% HPD

    outputContent += `clearance rate (IU/min)\n`;
    outputContent += `  Median: ${dd_median.toFixed(4)}\n`;
    outputContent += `  Mean: ${dd_mean.toFixed(4)}\n`;
    outputContent += `  95% HPD: [${dd_hpd_min.toFixed(4)}, ${dd_hpd_max.toFixed(4)}]\n\n`;
    outputContent += `frequency (pulses/h)\n`;
    outputContent += `  Median: ${freq_median.toFixed(4)}\n`;
    outputContent += `  Mean: ${freq_mean.toFixed(4)}\n`;
    outputContent += `  95% HPD: [${freq_hpd_min.toFixed(4)}, ${freq_hpd_max.toFixed(4)}]\n\n`;
  //  outputContent += `amplitude (IU)\n`;
  //  outputContent += `  Median: ${amplitude_median.toFixed(4)}\n`;
  //  outputContent += `  Mean: ${amplitude_mean.toFixed(4)}\n`;
  //  outputContent += `  95% HPD: [${amplitude_hpd_min.toFixed(4)}, ${amplitude_hpd_max.toFixed(4)}]\n\n`;
    
    
    // Write to file with the specified prefix
    //const outDir = path.join(outputfolder, 'out');
    const outputFileName = `${outputPrefix}_parameter_statistics.txt`;
    const outputPath = path.join(outputfolder, outputFileName);
    
    await fsp.writeFile(outputPath, outputContent);

    console.log(`Parameter statistics written to: ${outputPath}`);
    return outputPath;
  } catch (error) {
    console.error('Error in writeParamData:', error);
    throw error;
  }
}

function calculateMedian(data) {
  const sorted = [...data].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  if (sorted.length % 2 === 0) {
    return (sorted[middle - 1] + sorted[middle]) / 2;
  }
  return sorted[middle];
}

function calculateMean(data) {
  return data.reduce((sum, value) => sum + value, 0) / data.length;
}

// Make sure to add this IPC handler in your main process
ipcMain.handle('write-param-data', async (event, files, outputfolder, outputPrefix, burn_) => {
  try {
    const outputPath = await writeParamData(files, outputfolder, outputPrefix, burn_);
    return outputPath;
  } catch (error) {
    console.error('Error in write-param-data handler:', error);
    throw error;
  }
});

ipcMain.handle('write-file', async (event, fileName, outputfolder, content) => {
  
  //const outDir = path.join(outputfolder, 'out');
  const filePath = path.join(outputfolder, fileName);
  
  try {
    // Ensure the out directory exists
    //await fsp.mkdir(outputfolder, { recursive: true });
    
    // Write the file
    await fsp.writeFile(filePath, content);
    console.log(`File written successfully: ${filePath}`);
    return true;
  } catch (error) {
    console.error('Error writing file:', error);
    return false;
  }
});
/* HANDLING SETTINGS FILE  */

async function readSettingsFile() {
  const appPath = app.getAppPath();
  const settingsPath = path.join(appPath, 'data', 'params.txt');
  
  console.log('Attempting to read settings file:', settingsPath);
  
  try {
    const data = await fsp.readFile(settingsPath, 'utf8');
    console.log('Successfully read settings file. Content:', data);
    const lines = data.split('\n');
    const settings = {};

    lines.forEach((line, index) => {
      console.log(`Processing line ${index + 1}: ${line}`);
      const parts = line.trim().split(/\s+/);
      if (parts.length >= 4) {
        const [description, name, type, ...valueParts] = parts;
        const value = valueParts.join(' ');
        if (!settings[description]) {
          settings[description] = {};
        }
        settings[description][name] = { type, value };
        console.log(`Added setting: ${description}.${name} = ${value}`);
      } else {
        console.warn(`Skipping invalid line ${index + 1}: ${line}`);
      }
    });

    console.log('Parsed settings:', JSON.stringify(settings, null, 2));
    return settings;
  } catch (error) {
    console.error('Error reading settings file:', error);
    return null;
  }
}

ipcMain.handle('get-settings', async () => {
  console.log('Received get-settings request');
  const settings = await readSettingsFile();
  console.log('Returning settings:', settings);
  return settings;
});

ipcMain.handle('open-folder', async (event, folder) => {
  
  shell.openPath(folder);

});

ipcMain.handle('save-settings', async (event, newSettings) => {
  const appPath = app.getAppPath();
  const settingsPath = path.join(appPath, 'data', 'params.txt');
  
  console.log('Received save-settings request', newSettings);
  
  try {
    let content = '';
    for (const [group, settings] of Object.entries(newSettings)) {
      for (const [name, { type, value }] of Object.entries(settings)) {
        content += `${group} ${name} ${type} ${value}\n`;
      }
    }
    await fsp.writeFile(settingsPath, content.trim());
    
    logMessage(`Settings saved successfully in ${settingsPath}`);
    return true;
  } catch (error) {
    logMessage('An error occuered while saving settings');
    return false;
  }
});



/** plot mean LH and Gen data*/
ipcMain.handle('load-and-process-trajectory-data', async (event, files, burn_) => {
  const processedData = [];

  for (const file of files) {
    console.log('file', file);
    const fileContent = await fsp.readFile(file, 'utf8');
    const rows = fileContent.split('\n').map(row => row.split(','));

    const genData = [];
    const LHData = [];
    const timeData = [];

    const removeCount = Math.floor((rows.length/23) * (burn_ / 100));

    // Start from index 4 (5th row) for genData and index 6 (7th row) for LHData
    // Remember, array indexing starts at 0, so we subtract 1 from the row numbers
    let k =0;
    for (let i = 4; i < rows.length; i += 23) {
      if (rows[i] && k>=removeCount) genData.push(rows[i]);
      k = k+1;
    }

    k=0
    for (let i = 6; i < rows.length; i += 23) {
      if (rows[i] && k>=removeCount) LHData.push(rows[i]);
      k = k+1;
    }

    
    timeData.push(rows[22]);
    
    processedData.push({ timeData, genData, LHData });
    
  }

  
  
  return processedData;
});

ipcMain.handle('open-folder-dialog', async () => {
  const result = await dialog.showOpenDialog({
    properties: ['openDirectory']
  });
  return result;
});


ipcMain.handle('read-outputs', async (event, folder) => {
  try {
    return new Promise((resolve, reject) => {
      fs.readdir(folder, (err, files) => {
          if (err) {
              reject(err);
          } else {
              resolve(files);
          }
      });
    });
  } catch (error) {
    console.error('Error reading from out folder:', error);
    throw error;
  }
})

