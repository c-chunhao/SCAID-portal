import process from 'node:process'
import { fileURLToPath, URL } from 'node:url'

import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'

import AutoImport from 'unplugin-auto-import/vite'
import Components from 'unplugin-vue-components/vite'
import { ElementPlusResolver } from 'unplugin-vue-components/resolvers'

// https://vite.dev/config/
export default defineConfig({
  plugins: [
    vue(),
    AutoImport({
      resolvers: [ElementPlusResolver()],
    }),
    Components({
      resolvers: [
        ElementPlusResolver({ importStyle: "sass" })
      ],
    }),
  ],
  server: {
    host: '0.0.0.0',
    // Local development only: production serves /api and /system/media from the same Nginx origin.
    proxy: {
      '/api': {
        target: process.env.SCAID_API_TARGET || 'http://127.0.0.1:18080',
        changeOrigin: true,
      },
      '/system/media/': {
        target: process.env.SCAID_MEDIA_TARGET || 'http://127.0.0.1:10209',
        changeOrigin: true,
      },
    },
  },
  resolve: {
    alias: {
      '@': fileURLToPath(new URL('./src', import.meta.url))
    },
  },
  css: {
    preprocessorOptions: {
      scss: {
        additionalData: `
          @use "@/styles/element/index.scss" as *;
          @use "@/styles/var.scss" as *;
        `,
      },
    },
  },
})
