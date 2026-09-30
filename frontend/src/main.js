import '@/styles/bootstrap.scss'
import '@/styles/scaid.scss'

import { createApp } from 'vue'
import { createPinia } from 'pinia'

import App from './App.vue'
import router from './router'
import { lazyPlugin } from '@/directives/index';

import mitt from "mitt"


const app = createApp(App)

app.config.globalProperties.$bus = new mitt()

app.use(createPinia())
app.use(router)
app.use(lazyPlugin)

app.mount('#app')
