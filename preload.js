const { contextBridge, ipcRenderer } = require('electron');



ipcRenderer.on('app-version', (event, version) => {
  document.title = `HormoneBayes Analysis Tool (v${version})`;
});

//ipcRenderer.on('app-path', (event, appPath) => {
//  const viewResultsButton = document.getElementById('load-results');
//  viewResultsButton.setAttribute('data-tooltip', appPath);
//});
ipcRenderer.on('app-out-path', (event, appPath) => {
  
  const appFolderLabel = document.getElementById('app-folder-label');
  appFolderLabel.textContent = `${appPath}`;
  
});

ipcRenderer.on('app-data-path', (event, appPath) => {
  
  const appDataFolderLabel = document.getElementById('app-data-folder-label');
  appDataFolderLabel.textContent = `${appPath}`;

  const appSettingFolderLabel = document.getElementById('app-settings-folder-label');
  appSettingFolderLabel.textContent = `${appPath}`;


});

contextBridge.exposeInMainWorld('electron', {
  
  ipcRenderer: {
    send: (channel, ...args) => ipcRenderer.send(channel, ...args),
    invoke: (channel, ...args) => {
      console.log(`Invoking ${channel} with args:`, args);
      return ipcRenderer.invoke(channel, ...args);
    },
    on: (channel, func) => ipcRenderer.on(channel, (event, ...args) => func(event, ...args)),
    
  },
  
  openFolderDialog: () => ipcRenderer.invoke('open-folder-dialog'),
  getAppVersion: () => ipcRenderer.invoke('get-app-version')
});

window.addEventListener('DOMContentLoaded', () => {
  const replaceText = (selector, text) => {
    const element = document.getElementById(selector);
    if (element) element.innerText = text;
  };

  for (const dependency of ['chrome', 'node', 'electron']) {
    replaceText(`${dependency}-version`, process.versions[dependency]);
  }
});


